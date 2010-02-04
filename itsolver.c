#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include <mpi.h>

#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#else
#ifdef _FREEBSD
#include <cblas.h>
#endif
#endif

/* These headers are provided by ScaleME. */
#include "ScaleME.h"
#include "gmres.h"

#include "mlfma.h"
#include "direct.h"
#include "itsolver.h"

static int compcrt (complex float *dst, complex float *src) {
	int i;

#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < fmaconf.numbases; ++i)
		dst[i] = src[i] * fmaconf.contrast[i];

	return fmaconf.numbases;
}

static int augcrt (complex float *dst, complex float *src) {
	int i;

#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < fmaconf.numbases; ++i)
		dst[i] = src[i] - dst[i];

	return fmaconf.numbases;
}

int bicgstab (complex float *rhs, complex float *sol, int silent, solveparm *slv) {
	int i, j, rank;
	complex float *r, *rhat, *v, *p, *mvp, *t, *s;
	complex float rho, alpha, omega[2], beta, rnorm;
	float err, errinc;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	rho = alpha = omega[0] = 1.;

	/* Allocate the work arrays. */
	r = malloc (7 * fmaconf.numbases * sizeof(complex float));
	rhat = r + fmaconf.numbases;
	v = rhat + fmaconf.numbases;
	p = v + fmaconf.numbases;
	t = p + fmaconf.numbases;
	s = t + fmaconf.numbases;
	mvp = s + fmaconf.numbases;

	/* Zero the solution buffer and copy the residuals. */
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < fmaconf.numbases; ++i) {
		r[i] = rhat[i] = rhs[i];
		v[0] = sol[i] = 0;
	}

	/* Find the initial residual magnitude. */
	cblas_cdotc_sub (fmaconf.numbases, r, 1, r, 1, &rnorm);
	MPI_Allreduce(MPI_IN_PLACE, &rnorm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	errinc = sqrt(creal(rnorm));

	for (i = 0; i < slv->maxit; ++i) {
		/* Pre-compute portion of beta from previous iteration. */
		beta = alpha / (rho * omega[0]);
		/* Collect the total rho. */
		cblas_cdotc_sub (fmaconf.numbases, rhat, 1, r, 1, &rho);
		MPI_Allreduce(MPI_IN_PLACE, &rho, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		/* Include the missing factor. */
		beta *= rho;

#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < fmaconf.numbases; ++j) {
			/* Update the search vector. */
			p[j] = r[j] + beta * (p[j] - omega[0] * v[j]);
			/* Pre-compute the contrast for FMM evaluation. */
			mvp[j] = p[j] * fmaconf.contrast[j];
		}
		/* Matrix-vector product. */
		clrdircache ();
		ScaleME_applyParFMA(mvp, v, 0);
		/* Identity portion. */
		augcrt (v, p);
		/* Compute the inner product for alpha. */
		cblas_cdotc_sub (fmaconf.numbases, rhat, 1, v, 1, &alpha);
		MPI_Allreduce(MPI_IN_PLACE, &alpha, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		alpha = rho / alpha;
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < fmaconf.numbases; ++j) {
			s[j] = r[j] - alpha * v[j];
			/* Pre-compute the next contrast. */
			mvp[j] = s[j] * fmaconf.contrast[j];
		}
		/* Matrix-vector product. */
		clrdircache ();
		ScaleME_applyParFMA(mvp, t, 0);
		/* Identity portion. */
		augcrt (t, s);
		/* Compute the numerator of omega. */
		cblas_cdotc_sub (fmaconf.numbases, t, 1, s, 1, omega);
		/* Compute the denominator of omega. */
		cblas_cdotc_sub (fmaconf.numbases, t, 1, t, 1, omega + 1);
		MPI_Allreduce(MPI_IN_PLACE, omega, 4, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		/* Compute the ratio for omega. */
		omega[0] /= omega[1];
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < fmaconf.numbases; ++j) {
			/* Update the solution vector. */
			sol[j] += alpha * p[j] + omega[0] * s[j];
			/* Update the residual vector. */
			r[j] = s[j] - omega[0] * t[j];
		}

		/* Compute the norm of the residual s. */
		cblas_cdotc_sub (fmaconf.numbases, r, 1, r, 1, &rnorm);
		MPI_Allreduce(MPI_IN_PLACE, &rnorm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		err = sqrt(creal(rnorm)) / errinc;

		if (!rank && !silent) printf ("BiCG-STAB(%d): %g\n", i, err);

		/* Break if convergence is detected. */
		if (err < slv->epscg) break;
	}

	free (r);

	return i;
}

int cgmres (complex float *rhs, complex float *sol, int silent, solveparm *slv) {
	int icntl[8], irc[5], lwork, info[3], i, myRank;
	float rinfo[2], cntl[5];
	complex float *zwork, *solbuf, *tx, *ty, *tz;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	/* Allocate memory for work array. */
	lwork = slv->restart * slv->restart +
		slv->restart * (fmaconf.numbases + 5) + 5 * fmaconf.numbases + 2;
	zwork = calloc (lwork, sizeof(complex float));
	solbuf = malloc (fmaconf.numbases * sizeof(complex float));

	/* Initialize the parameters. */
	initcgmres_(icntl, cntl);

	/* Only the root process should print convergence history. */
	if (!myRank && !silent) icntl[0] = icntl[1] = icntl[2] = 6;
	else icntl[1] = icntl[2] = 0;

	/* None of the processes should make noise about GMRES. */
	if (silent) icntl[0] = icntl[1] = icntl[2] = 0;

	/* A preconditioner will not be used. */
	icntl[3] = 0;

	icntl[4] = 0; /* Use MGS for orthogonalization. */
	icntl[5] = 1; /* Use an initial guess: the incident field. */
	icntl[6] = slv->maxit; /* Set the maximum interation count. */

	cntl[0] = slv->epscg;

	/* Copy the initial guess: use the RHS. */
	memcpy (zwork, rhs, fmaconf.numbases * sizeof(complex float));
	/* Copy the RHS. */
	memcpy (zwork + fmaconf.numbases, rhs, fmaconf.numbases * sizeof(complex float));

	do {
		drivecgmres_(&(fmaconf.gnumbases), &(fmaconf.numbases),
				&(slv->restart), &lwork, zwork, irc,
				icntl, cntl, info, rinfo);
		if (!(info[0]) && !(irc[0])) break;

		switch (irc[0]) {
		case EXIT: case RIGHT_PRECOND: case LEFT_PRECOND: break;
		case MATVEC:
			   /* fortran indices start from 1 */
			   tx = zwork+irc[1] - 1;
			   ty = zwork+irc[3] - 1;
			   clrdircache ();
			   compcrt (solbuf, tx);
			   ScaleME_applyParFMA(solbuf, ty, 0);
			   augcrt (ty, tx);
			   break;
		case DOT_PROD:
			   ty = zwork + irc[2] - 1;
			   tx = zwork + irc[1] - 1;
			   tz = zwork + irc[3] - 1;

#pragma omp parallel for default(shared) private(i)
			   for (i = 0; i < irc[4]; ++i)
				   /* compute the local dot product first */
				   cblas_cdotc_sub (fmaconf.numbases, tx + i * fmaconf.numbases, 1, ty, 1, tz + i);
			   
			   /* now do a global reduce to get the final answer */
			   MPI_Allreduce(MPI_IN_PLACE, tz, 2 * irc[4], MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
			   break;
		}
	} while (irc[0]);

	if (!myRank && info[0]) fprintf (stdout, "CGMRES: return value: %d\n", info[0]); 

	memcpy (sol, zwork, fmaconf.numbases * sizeof(complex float));
	
	if (!myRank && !silent)
		fprintf(stdout, "CGMRES: %d iterations, %.6E PBE, %.6E BE.\n", info[1], rinfo[0], rinfo[1]);

	free (zwork);
	free (solbuf);
	return info[0];
}

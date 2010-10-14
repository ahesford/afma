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

int matvec (complex float *out, complex float *in, complex float *cur) {
	long i, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;

	/* Compute the contrast pressure. */
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < nelt; ++i) cur[i] = in[i] * fmaconf.contrast[i];

	/* Reset the direct-interaction buffer and compute
	 * the matrix-vector product for the Green's matrix. */
	clrdircache();
	ScaleME_applyParFMA (cur, out);

	/* Add in the identity portion. */
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < nelt; ++i) out[i] = in[i] - out[i];

	return 0;
}

static complex float pardot (complex float *x, complex float *y, long n) {
	complex float dp;

	/* Compute the local portion. */
	cblas_cdotc_sub (n, x, 1, y, 1, &dp);
	/* Add in the contributions from other processors. */
	MPI_Allreduce (MPI_IN_PLACE, &dp, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return dp;
}

int bicgstab (complex float *rhs, complex float *sol,
		int guess, int mit, float tol, int quiet) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	int i, rank;
	complex float *r, *rhat, *v, *p, *mvp, *t;
	complex float rho, alpha, omega, beta;
	float err, rhn;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	rho = alpha = omega = 1.;

	/* Allocate and zero the work arrays. */
	r = calloc (6L * nelt, sizeof(complex float));
	rhat = r + nelt;
	v = rhat + nelt;
	p = v + nelt;
	t = p + nelt;
	mvp = t + nelt;

	/* Compute the norm of the right-hand side for residual scaling. */
	rhn = sqrt(creal(pardot (rhs, rhs, nelt)));

	/* Compute the inital matrix-vector product for the input guess. */
	if (guess) matvec (r, sol, mvp);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < nelt; ++j) r[j] = rhs[j] - r[j];

	if (!guess) memset (sol, 0, nelt * sizeof(complex float));
		
	/* Copy the initial residual as the test vector. */
	memcpy (rhat, r, nelt * sizeof(complex float));

	/* Find the norm of the initial residual. */
	err = sqrt(creal(pardot (r, r, nelt))) / rhn;
	if (!rank && !quiet) printf ("True residual: %g\n", err);

	/* Run iterations until convergence or the maximum is reached. */
	for (i = 0; i < mit && err > tol; ++i) {
		/* Pre-compute portion of beta from previous iteration. */
		beta = alpha / (rho * omega);
		/* Compute rho for this iteration. */
		rho = pardot (rhat, r, nelt);
		/* Include the missing factor in beta. */
		beta *= rho;

		/* Update the search vector. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j)
			p[j] = r[j] + beta * (p[j] - omega * v[j]);

		/* Compute the first search step, v = A * p. */
		matvec (v, p, mvp);

		/* Compute the next alpha. */
		alpha = rho / pardot (rhat, v, nelt);

#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) {
			/* Update the solution vector. */
			sol[j] += alpha * p[j];
			/* Update the residual vector. */
			r[j] -= alpha * v[j];
		}

		/* Compute the scaled residual norm and stop if convergence
		 * has been achieved. */
		err = sqrt(creal(pardot (r, r, nelt))) / rhn;
		if (!rank && !quiet) printf ("BiCG-STAB(%0.1f): %g\n", 0.5 + i, err);

		/* Flush the output buffers. */
		fflush (stdout);
		fflush (stderr);

		if (err < tol) break;

		/* Compute the next search step, t = A * r. */
		matvec (t, r, mvp);

		/* Compute the update direction. */
		omega = pardot (t, r, nelt) / pardot (t, t, nelt);

		/* Update both the residual and the solution guess. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) {
			/* Update the solution vector. */
			sol[j] += omega * r[j];
			/* Update the residual vector. */
			r[j] -= omega * t[j];
		}
	
		/* Compute the scaled residual norm. */
		err = sqrt(creal(pardot (r, r, nelt))) / rhn;
		if (!rank && !quiet) printf ("BiCG-STAB(%d): %g\n", i + 1, err);
	}

	free (r);
	return i;
}

int cgmres (complex float *rhs, complex float *sol,
		int guess, int silent, solveparm *slv) {
	int icntl[8], irc[5], lwork, info[3], i, myRank;
	int nelt = fmaconf.numbases * fmaconf.bspboxvol;
	float rinfo[2], cntl[5];
	complex float *zwork, *solbuf, *tx, *ty, *tz;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	/* Allocate memory for work array. */
	lwork = slv->restart * slv->restart +
		slv->restart * (nelt + 5) + 5 * nelt + 2;
	zwork = calloc (lwork, sizeof(complex float));
	solbuf = malloc (nelt * sizeof(complex float));

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
	icntl[5] = guess; /* Use an initial guess: the incident field. */
	icntl[6] = slv->maxit; /* Set the maximum interation count. */

	cntl[0] = slv->epscg;

	/* Copy the initial guess, if it exists; otherwise zero the workspace. */
	if (guess) memcpy (zwork, rhs, nelt * sizeof(complex float));
	else memset (zwork, 0, nelt * sizeof(complex float));

	/* Copy the RHS. */
	memcpy (zwork + nelt, rhs, nelt * sizeof(complex float));

	do {
		drivecgmres_(&(fmaconf.gnumbases), &(nelt),
				&(slv->restart), &lwork, zwork, irc,
				icntl, cntl, info, rinfo);
		if (!(info[0]) && !(irc[0])) break;

		switch (irc[0]) {
		case EXIT: case RIGHT_PRECOND: case LEFT_PRECOND: break;
		case MATVEC:
			   /* fortran indices start from 1 */
			   tx = zwork+irc[1] - 1;
			   ty = zwork+irc[3] - 1;
			   matvec (ty, tx, solbuf);
			   break;
		case DOT_PROD:
			   ty = zwork + irc[2] - 1;
			   tx = zwork + irc[1] - 1;
			   tz = zwork + irc[3] - 1;

#pragma omp parallel for default(shared) private(i)
			   for (i = 0; i < irc[4]; ++i)
				   /* compute the local dot product first */
				   cblas_cdotc_sub (nelt, tx + i * nelt, 1, ty, 1, tz + i);
			   
			   /* now do a global reduce to get the final answer */
			   MPI_Allreduce(MPI_IN_PLACE, tz, 2 * irc[4], MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
			   break;
		}

		/* Flush output buffers to report residual errors. */
		fflush (stdout);
		fflush (stderr);
	} while (irc[0]);

	if (!myRank && info[0]) printf ("CGMRES: return value: %d\n", info[0]); 

	memcpy (sol, zwork, nelt * sizeof(complex float));
	
	if (!myRank && !silent)
		printf("CGMRES: %d iterations, %.6E PBE, %.6E BE.\n",
				info[1], rinfo[0], rinfo[1]);

	free (zwork);
	free (solbuf);
	return info[1];
}

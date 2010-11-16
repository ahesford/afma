#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include <mpi.h>

/* Pull in the CBLAS header. */
#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#else
#ifdef _FREEBSD
#include <cblas.h>
#else
#include <gsl_cblas.h>
#endif /* _FREEBSD */
#endif /* _MACOSX */

/* These headers are provided by ScaleME. */
#include "ScaleME.h"

#include "mlfma.h"
#include "direct.h"
#include "itsolver.h"
#include "util.h"

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

int gmres (complex float *rhs, complex float *sol,
		int guess, int mit, float tol, int quiet) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol, lwork;
	int i, rank, one = 1;
	complex float *h, *v, *mvp, *beta, *y;
	complex float *vp, *hp, *s, cr, cone = 1.;
	float rhn, err, *c;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Allocate space for all required complex vectors. */
	lwork = (mit + 1) * (mit + nelt + 1) + nelt + mit;
	v = calloc (lwork, sizeof(complex float));	/* The Krylov subspace. */
	beta = v + nelt * (mit + 1);		/* The least-squares RHS. */
	mvp = beta + mit + 1;			/* Buffer for matrix-vector product. */
	h = mvp + nelt;				/* The upper Hessenberg matrix. */
	s = h + (mit + 1) * mit;		/* Givens rotation sines. */

	/* Allocate space for the Givens rotation cosines. */
	c = malloc (mit * sizeof(float));

	/* Compute the norm of the RHS for residual scaling. */
	rhn = parnorm(rhs, nelt);

	/* Compute the initial matrix-vector product for the input guess. */
	if (guess) matvec (v, sol, mvp);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < nelt; ++j) v[j] = rhs[j] - v[j];

	/* Zero the initial guess if one wasn't provided. */
	if (!guess) memset (sol, 0, nelt * sizeof(complex float));

	/* Find the norm of the initial residual. */
	err = parnorm(v, nelt);

	/* Construct the initial Arnoldi vector by normalizing the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < nelt; ++j) v[j] /= err;

	/* Construct the vector beta for the minimization problem. */
	beta[0] = err;

	/* Report the RRE. */
	err /= rhn;
	if (!rank && !quiet) printf ("True residual: %g\n", err);

	for (i = 0; i < mit && err > tol; ++i) {
		/* Point to the working space for this iteration. */
		vp = v + i * nelt;
		hp = h + i * (mit + 1);

		/* Compute the next expansion of the Krylov space. */
		matvec (vp + nelt, vp, mvp);
		/* Perform modified Gram-Schmidt to orthogonalize the basis. */
		/* This also builds the Hessenberg matrix column. */
		cmgs (vp + nelt, hp, v, nelt, i + 1);
		/* Compute the norm of the next basis vector. */
		hp[i + 1] = parnorm(vp + nelt, nelt);

		/* Avoid breakdown. */
		if (cabs(hp[i + 1]) <  FLT_EPSILON) break;

		/* Normalize the basis vector. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) vp[nelt + j] /= creal(hp[i + 1]);

		/* Apply previous Givens rotations to the Hessenberg column. */
		for (j = 0; j < i; ++j) 
			crot_ (&one, hp + j, &one, hp + j + 1, &one, c + j, s + j);

		/* Compute the Givens rotation for the current iteration. */
		clartg_ (hp + i, hp + i + 1, c + i, s + i, &cr);
		/* Apply the current Givens rotation to the Hessenberg column. */
		hp[i] = cr;
		hp[i + 1] = 0;
		/* Perform the rotation on the vector beta. */
		crot_ (&one, beta + i, &one, beta + i + 1, &one, c + i, s + i);

		/* Estimate the RRE for this iteration. */
		err = cabs(beta[i + 1]) / rhn;
		if (!rank && !quiet) printf ("GMRES(%d): %g\n", i, err);
	}

	/* If there were any GMRES iterations, update the solution. */
	if (i > 0) {
		/* Compute the minimizer of the least-squares problem. */
		cblas_ctrsv (CblasColMajor, CblasUpper, CblasNoTrans,
				CblasNonUnit, i, h, mit + 1, beta, 1);
		
		/* Compute the update to the solution. */
		cblas_cgemv (CblasColMajor, CblasNoTrans, nelt, i,
				&cone, v, nelt, beta, 1, &cone, sol, 1);
	}

	free (v);
	free (c);

	return i;
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
	rhn = parnorm(rhs, nelt);

	/* Compute the inital matrix-vector product for the input guess. */
	if (guess) matvec (r, sol, mvp);

	/* Subtract from the RHS to form the residual. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < nelt; ++j) r[j] = rhs[j] - r[j];

	if (!guess) memset (sol, 0, nelt * sizeof(complex float));
		
	/* Copy the initial residual as the test vector. */
	memcpy (rhat, r, nelt * sizeof(complex float));

	/* Find the norm of the initial residual. */
	err = parnorm(r, nelt) / rhn;
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
		err = parnorm(r, nelt) / rhn;
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
		err = parnorm(r, nelt) / rhn;
		if (!rank && !quiet) printf ("BiCG-STAB(%d): %g\n", i + 1, err);

		fflush (stdout);
		fflush (stderr);
	}

	free (r);
	return i;
}

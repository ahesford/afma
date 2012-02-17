#include <string.h>

#include <mpi.h>

#include "ScaleME.h"

#include "precision.h"

#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"
#include "frechet.h"

/* Computes a contribution to the Frechet derivative for one transmitter. */
int frechet (cplx *crt, cplx *fld,
		cplx *sol, measdesc *obs, solveparm *slv) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	cplx *zwork, *zwcrt;

	/* Set up the workspaces. */
	zwork = malloc (2L * nelt * sizeof(cplx));
	zwcrt = zwork + nelt;

	/* Given the test vector and field distribution, compute the currents. */
	for (j = 0; j < nelt; ++j) zwcrt[j] = crt[j] * fld[j];

	/* Compute the RHS for the given current distribution.
	 * Store it in the buffer space. */
	ScaleME_applyParFMA (zwcrt, zwork);

	/* Compute the Frechet derivative field. */
	bicgstab (zwork, zwork, 0, slv->maxit, slv->epscg, 1);

	/* Convert this into a current distribution. */
	for (j = 0; j < nelt; ++j)
		zwork[j] = fmaconf.contrast[j] * zwork[j] + zwcrt[j];

	/* Compute the measured scattered field. */
	farfield (zwork, obs, sol);

	free (zwork);
	return obs->count;
}

int frechadj (cplx *mag, cplx *fld,
		cplx *sol, measdesc *obs, solveparm *slv) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	real factor = fmaconf.k0 * fmaconf.k0 * fmaconf.cellvol;
	cplx *zwork, *smag, scale;

	/* Temporary workspaces. */
	zwork = malloc ((nelt + (long)obs->count) * sizeof(cplx));
	smag = zwork + nelt;

	/* Compensate incident field for FMM incoming signature scaling. */
	scale = I * fmaconf.k0 * fmaconf.k0;

	/* Conjugate for back-propagtion. Scale to compensate for
	 * scaling of incoming far-field signature. */
	for (j = 0; j < obs->count; ++j)
		smag[j] = conj(mag[j]) / scale;

	/* Compute the RHS for the provided magnitude distribution. */
	ScaleME_setRootFarFld (obs->imat[1], zwork, &smag);

	/* Compute the adjoint Frechet derivative field. */
	bicgstab (zwork, zwork, 0, slv->maxit, slv->epscg, 1);

	/* Augment the solution for this transmitter. */
	for (j = 0; j < nelt; ++j)
		sol[j] += factor * conj (zwork[j] * fld[j]);

	free (zwork);

	return fmaconf.numbases;
}

/* Estimate the square of the spectral radius of the Frechet derivative
 * using power iteration on the product of the derivative and its adjoint. */
real specrad (int maxit, solveparm *slv, measdesc *src, measdesc *obs) {
	cplx *ifld, *adjcrt, *scat, *pn, mu = 1.0, newmu;
	long k, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol,
	     randmax = (2L << 31) - 1;
	real nrm, tmp;
	int i, j;

	ifld = malloc (3L * nelt * sizeof(cplx));
	adjcrt = ifld + nelt;
	pn = adjcrt + nelt;

	scat = malloc (obs->count * sizeof(cplx));

	/* Build a random initial vector. */
	for (k = 0, nrm = 0.0; k < nelt; ++k) {
		pn[k] = ((real)random() + I * (real)random()) / (real)randmax;
		tmp = cabs(pn[k]);
		nrm += tmp * tmp;
	}
	MPI_Allreduce (MPI_IN_PLACE, &nrm, 1, MPIREAL, MPI_SUM, MPI_COMM_WORLD);
	nrm = sqrt(nrm);

	/* Normalize the initial vector. */
	for (k = 0; k < nelt; ++k) pn[k] /= nrm;

	/* Loop through power iterations. */
	for (j = 0; j < maxit; ++j) {
		memset (adjcrt, 0, nelt * sizeof(cplx));

		/* Compute the product (F*)(F)p. */
		for (i = 0; i < src->count; ++i) {
			/* Build the incident field. */
			buildrhs (ifld, src->locations + 3 * i, src->plane, NULL);
			/* Solve for the internal field. */
			bicgstab (ifld, ifld, 0, slv->maxit, slv->epscg, 1);
			/* Compute the Frechet derivative. */
			frechet (pn, ifld, scat, obs, slv);
			/* Compute the adjoint Frechet derivative. */
			frechadj (scat, ifld, adjcrt, obs, slv);
		}

		/* Update the norm of the solution and find the approximate
		 * singular value. Note that the norm of the test vector is 1. */
		for (k = 0, nrm = newmu = 0.0; k < nelt; ++k) {
			newmu += conj(pn[k]) * adjcrt[k];
			tmp = cabs(adjcrt[k]);
			nrm += tmp * tmp;
		}
		MPI_Allreduce (MPI_IN_PLACE, &nrm, 1, MPIREAL, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce (MPI_IN_PLACE, &newmu, 2, MPIREAL, MPI_SUM, MPI_COMM_WORLD);
		nrm = sqrt(nrm);

		/* Break if the updated singular value hasn't changed much. */
		if (cabs((newmu - mu) / mu) < slv->epscg) break;

		/* Normalize the solution for the next iteration. */
		for (k = 0; k < nelt; ++k)
			pn[k] = adjcrt[k] / nrm;

		/* Copy the new singular value. */
		mu = newmu;
	}

	free (ifld);
	free (scat);

	return creal(mu);
}

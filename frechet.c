#include <math.h>
#include <complex.h>

#include <mpi.h>

#include "ScaleME.h"

#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"
#include "frechet.h"

/* Computes a contribution to the Frechet derivative for one transmitter. */
int frechet (complex float *crt, complex float *fld,
		complex float *sol, measdesc *obs, solveparm *slv) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	complex float *zwork, *zwcrt;

	/* Set up the workspaces. */
	zwork = malloc (2L * nelt * sizeof(complex float));
	zwcrt = zwork + nelt;

	/* Given the test vector and field distribution, compute the currents. */
	for (j = 0; j < nelt; ++j) zwcrt[j] = crt[j] * fld[j];

	/* Compute the RHS for the given current distribution.
	 * Store it in the buffer space. */
	ScaleME_applyParFMA (zwcrt, zwork);

	/* Compute the Frechet derivative field. */
	cgmres (zwork, zwork, 1, slv);

	/* Convert this into a current distribution. */
	for (j = 0; j < nelt; ++j)
		zwork[j] = fmaconf.contrast[j] * zwork[j] + zwcrt[j];

	/* Compute the measured scattered field. */
	farfield (zwork, obs, sol);

	free (zwork);
	return obs->count;
}

int frechadj (complex float *mag, complex float *fld,
		complex float *sol, measdesc *obs, solveparm *slv) {
	long j, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	float factor = fmaconf.k0 * fmaconf.k0 * fmaconf.cellvol;
	complex float *zwork, *smag, scale;

	/* Temporary workspaces. */
	zwork = malloc ((nelt + (long)obs->count) * sizeof(complex float));
	smag = zwork + nelt;

	/* Compensate incident field for FMM incoming signature scaling. */
	scale = I * fmaconf.k0 * fmaconf.k0;

	/* Conjugate for back-propagtion. Scale to compensate for
	 * scaling of incoming far-field signature. */
	for (j = 0; j < obs->count; ++j)
		smag[j] = conj(mag[j]) / scale;

	/* Compute the RHS for the provided magnitude distribution. */
	ScaleME_setRootFarFld (obs->imat, zwork, &smag);

	/* Compute the adjoint Frechet derivative field. */
	cgmres (zwork, zwork, 1, slv);

	/* Augment the solution for this transmitter. */
	for (j = 0; j < nelt; ++j)
		sol[j] += factor * conj (zwork[j] * fld[j]);

	free (zwork);

	return fmaconf.numbases;
}

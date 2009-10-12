#include <math.h>
#include <complex.h>

#include <mpi.h>

#include "ScaleME.h"

#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"
#include "frechet.h"
#include "excite.h"

complex float *zwork, *zwcrt, *rcvgrf;

complex float *bldfrechbuf (int size, measdesc *obs) {
	zwork = malloc (2 * size * sizeof(complex float));
	zwcrt = zwork + size;

	rcvgrf = malloc (fmaconf.numbases * obs->count * sizeof(complex float));
	precompgrf (obs, rcvgrf);

	return zwork;
}

void delfrechbuf (void) {
	if (zwork) free (zwork);
	if (rcvgrf) free (rcvgrf);
}

/* Computes a contribution to the Frechet derivative for one transmitter. */
int frechet (complex float *crt, complex float *fld,
		complex float *sol, measdesc *obs, solveparm *slv) {
	int j;

	/* Given the test vector and field distribution, compute the currents. */
	for (j = 0; j < fmaconf.numbases; ++j)
		zwcrt[j] = crt[j] * fld[j];

	/* Compute the RHS for the given current distribution.
	 * Store it in the buffer space. */
	ScaleME_applyParFMA (zwcrt, zwork, 0);

	/* Compute the Frechet derivative field. */
	cgmres (zwork, zwork, 1, slv);

	/* Convert this into a current distribution. */
	for (j = 0; j < fmaconf.numbases; ++j)
		zwork[j] = fmaconf.contrast[j] * zwork[j] + zwcrt[j];

	/* Compute the measured scattered field. */
	farfield (zwork, obs, sol);

	return obs->count;
}

int frechadj (complex float *mag, complex float *fld,
		complex float *sol, measdesc *obs, solveparm *slv) {
	int j;
	double factor = fmaconf.k0 * fmaconf.k0
		* fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];

	for (j = 0; j < obs->count; ++j)
		mag[j] = conj(mag[j]);

	/* Compute the RHS for the provided magnitude distribution.
	 * For now, the slow, direct calculation routine will be used. */
	precomprhs (zwork, obs, mag, rcvgrf);

	/* Compute the adjoint Frechet derivative field. */
	cgmres (zwork, zwork, 1, slv);

	/* Augment the solution for this transmitter. */
	for (j = 0; j < fmaconf.numbases; ++j)
		sol[j] += factor * conj (zwork[j] * fld[j]);

	for (j = 0; j < obs->count; ++j)
		mag[j] = conj(mag[j]);

	return fmaconf.numbases;
}

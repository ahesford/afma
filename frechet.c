#include <complex.h>

#include <mpi.h>

#include <Complex.h>
#include <ScaleME.h>

#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"
#include "frechet.h"
#include "excite.h"

complex float *zwork, *zwcrt;

complex float *bldfrechbuf (int size) {
	zwork = malloc (2 * size * sizeof(complex float));
	zwcrt = zwork + size;

	return zwork;
}

void delfrechbuf (void) {
	if (zwork) free (zwork);
}

/* Computes a contribution to the Frechet derivative for one transmitter. */
int frechet (complex float *crt, complex float *fld, complex float *sol) {
	int j;

	/* Given the test vector and field distribution, compute the currents. */
	for (j = 0; j < fmaconf.numbases; ++j)
		zwcrt[j] = crt[j] * fld[j];

	/* Compute the RHS for the given current distribution.
	 * Store it in the buffer space. */
	ScaleME_applyParFMA (REGULAR, zwcrt, zwork);

	/* Compute the Frechet derivative field. */
	cgmres (zwork, zwork);

	/* Convert this into a current distribution. */
	for (j = 0; j < fmaconf.numbases; ++j)
		zwork[j] = fmaconf.contrast[j] * zwork[j] + zwcrt[j];

	/* Compute the measured scattered field. */
	farfield (zwork, &obsmeas, sol);
	MPI_Bcast (sol, 2 * obsmeas.count, MPI_FLOAT, 0, MPI_COMM_WORLD);

	return obsmeas.count;
}

int frechadj (complex float *mag, complex float *fld, complex float *sol) {
	int j;

	for (j = 0; j < obsmeas.count; ++j)
		mag[j] = conj(mag[j]);

	/* Compute the RHS for the provided magnitude distribution.
	 * For now, the slow, direct calculation routine will be used. */
	multirhs (zwork, &obsmeas, mag, 1); 
	
	/* Compute the adjoint Frechet derivative field. */
	cgmres (zwork, zwork);

	/* Augment the solution for this transmitter. */
	for (j = 0; j < fmaconf.numbases; ++j)
		sol[j] += conj (zwork[j] * fld[j]);

	for (j = 0; j < obsmeas.count; ++j)
		mag[j] = conj(mag[j]);

	return fmaconf.numbases;
}

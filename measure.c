#include <string.h>
#include <complex.h>

#include <mpi.h>

#include <Complex.h> // ScaleME-provided complex include. */
#include <ScaleME.h>

#include "measure.h"
#include "mlfma.h"
#include "fsgreen.h"

measdesc srcmeas, obsmeas;

int farfield (complex float *currents, int nthetas, int nphis,
		float *thetas, float *phis, complex float *result) {
	Complex *fields = (Complex *)result;
	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, nthetas, nphis, thetas, phis,
				(Complex *)currents, &fields)) return 0;

	return 1;
}

int directfield (complex float *currents, int numobs,
		float *locations, complex float *result) {
	int i, j, gi;
	complex float val, *buf;
	float *cen;

	buf = malloc (numobs * sizeof(complex float));
	memset (buf, 0, numobs * sizeof(complex float));

	/* Compute the fields radiated by the local fields. */
	for (j = 0; j < fmaconf.numbases; ++j) {
		gi = fmaconf.bslist[j];
		bscenter (gi, cen);
		for (i = 0; i < numobs; ++i) {
			val = fsgreen (fmaconf.k0, cen, locations + 3 * i);
			buf[i] += val * currents[j];
		}
	}

	/* Multiply by cell volume. */
	for (i = 0; i < numobs; ++i)
		buf[i] *= fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];

	/* Collect the solution across all processors. */
	MPI_Allreduce (result, buf, 2 * numobs, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return numobs;
}

int buildlocs (measdesc *desc) {
	int i, j, k;
	float dtheta, dphi, theta, phi, rst;

	desc->count = desc->ntheta * desc->nphi;

	dtheta = (desc->trange[1] - desc->trange[0]);
	dtheta /= max (desc->ntheta - 1, 1);

	dphi = (desc->prange[1] - desc->prange[0]);
	dphi /= max (desc->nphi - 1, 1);

	desc->locations = malloc (3 * desc->count * sizeof(float));

	theta = desc->trange[0];
	for (i = 0, k = 0; i < desc->ntheta; ++i) {
		phi = desc->prange[0];
		for (j = 0; j < desc->nphi; ++j, ++k) {
			rst = desc->radius * sin (theta);
			desc->locations[3 * k] = rst * cos (phi);
			desc->locations[3 * k + 1] = rst * sin (phi);
			desc->locations[3 * k + 2] = desc->radius * cos (theta);
			phi += dphi;
		}
		theta += dtheta;
	}

	return desc->count;
}

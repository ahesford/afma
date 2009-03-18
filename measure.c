#include <string.h>
#include <math.h>
#include <complex.h>

#include <mpi.h>

#include <Complex.h> // ScaleME-provided complex include. */
#include <ScaleME.h>

#include "measure.h"
#include "mlfma.h"
#include "fsgreen.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

measdesc srcmeas, obsmeas;

int farfield (complex float *currents, measdesc *obs, complex float *result) {
	int i, j;
	float *thetas, dtheta;
	Complex *fields = (Complex *)result;
	complex float fact;

	thetas = malloc (obs->ntheta * sizeof(float));

	dtheta = (obs->trange[1] - obs->trange[0]);
	dtheta /= MAX(obs->ntheta - 1, 1);

	for (i = 0; i < obs->ntheta; ++i)
		thetas[i] = obs->trange[0] + i * dtheta;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, obs->ntheta, obs->nphi, thetas,
				obs->prange, (Complex *)currents, &fields))
		return 0;

	/* The far-field pattern already has a factor of k in the front.
	 * However, the actual integral needs (k^2 / 4 pi), so we need the
	 * extra factors in the field. */
	j = obs->ntheta * obs->nphi;

	/* Scale the pattern to produce actual measurements on the desired sphere. */
	fact = fmaconf.k0 * cexp (I * fmaconf.k0 * obs->radius) / (4 * M_PI * obs->radius);
	for (i = 0; i < j; ++i) result[i] *= fact;

	free (thetas);
	return 1;
}

int directfield (complex float *currents, measdesc *obs, complex float *result) {
	int i, j, gi;
	complex float val, *buf;
	float cen[3];

	buf = malloc (obs->count * sizeof(complex float));
	memset (buf, 0, obs->count * sizeof(complex float));

	/* Compute the fields radiated by the local fields. */
	for (j = 0; j < fmaconf.numbases; ++j) {
		gi = fmaconf.bslist[j];
		bscenter (gi, cen);
		for (i = 0; i < obs->count; ++i) {
			val = fsgreen (fmaconf.k0, cen, obs->locations + 3 * i);
			buf[i] += val * currents[j];
		}
	}

	for (i = 0; i < obs->count; ++i) {
		/* Scale by cell volume for one-point integration. */
		buf[i] *= fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];
		/* Include the factor of k0^2 in front of the integral. */
		buf[i] *= fmaconf.k0 * fmaconf.k0;
	}

	/* Collect the solution across all processors. */
	MPI_Allreduce (buf, result, 2 * obs->count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return obs->count;
}

int buildlocs (measdesc *desc) {
	int i, j, k;
	float dtheta, dphi, theta, phi, rst;

	desc->count = desc->ntheta * desc->nphi;

	dtheta = (desc->trange[1] - desc->trange[0]);
	dtheta /= MAX (desc->ntheta - 1, 1);

	dphi = (desc->prange[1] - desc->prange[0]);
	dphi /= MAX (desc->nphi - 1, 1);

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

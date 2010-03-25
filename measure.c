#include <string.h>
#include <math.h>
#include <complex.h>

#include <mpi.h>

#include "ScaleME.h"

#include "measure.h"
#include "mlfma.h"
#include "fsgreen.h"

#include "util.h"

int farfield (complex float *currents, measdesc *obs, complex float *result) {
	int i, mpirank;
	float *thetas, dtheta;
	complex float fact;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* The theta samples. */
	thetas = malloc (obs->ntheta * sizeof(float));

	dtheta = (obs->trange[1] - obs->trange[0]);
	dtheta /= MAX(obs->ntheta - 1, 1);

	for (i = 0; i < obs->ntheta; ++i)
		thetas[i] = obs->trange[0] + i * dtheta;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, obs->ntheta, obs->nphi, thetas,
				obs->prange, currents, &result)) return 0;

	/* The far-field pattern already has a factor of k in the front.
	 * However, the actual integral needs (k^2 / 4 pi), so we need the
	 * extra factors in the field. */
	fact = fmaconf.k0 / (4 * M_PI);
	for (i = 0; i < obs->count; ++i) result[i] *= fact;

	free (thetas);
	return 1;
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

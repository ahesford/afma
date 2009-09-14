#include <string.h>
#include <math.h>
#include <complex.h>

#include <mpi.h>

#include "ScaleME.h"
#include "ScaleME_Complex.h"

#include "measure.h"
#include "mlfma.h"
#include "fsgreen.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

int farfield (complex float *currents, measdesc *obs, complex float *result) {
	int i, mpirank;
	float *thetas, dtheta;
	complex float *fields = NULL, fact;
	Complex *fieldstct;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* The theta samples. */
	thetas = malloc (obs->ntheta * sizeof(float));

	if (!mpirank) fields = malloc (obs->count * sizeof(complex float));
	fieldstct = (Complex *)fields;

	dtheta = (obs->trange[1] - obs->trange[0]);
	dtheta /= MAX(obs->ntheta - 1, 1);

	for (i = 0; i < obs->ntheta; ++i)
		thetas[i] = obs->trange[0] + i * dtheta;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, obs->ntheta, obs->nphi, thetas,
				obs->prange, (Complex *)currents, &fieldstct))
		return 0;

	if (!mpirank) {
		/* The far-field pattern already has a factor of k in the
		 * front. However, the actual integral needs (k^2 / 4 pi), so
		 * we need the extra factors in the field. */
		fact = fmaconf.k0 / (4 * M_PI);
		for (i = 0; i < obs->count; ++i) result[i] = fact * fields[i];
	}

	/* Distribute the solution to all processors. */
	MPI_Bcast (result, 2 * obs->count, MPI_FLOAT, 0, MPI_COMM_WORLD);

	free (thetas);
	if (!mpirank) free (fields);
	return 1;
}

int directfield (complex float *currents, measdesc *obs,
		complex float *result, complex float *grf) {
	int i, j;
	complex float val, *buf, *gptr;
	float cen[3], fact;

	buf = calloc (obs->count, sizeof(complex float));

	/* The integration factor. */
	fact = fmaconf.k0 * fmaconf.k0 
		* fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];

	if (!grf) {
		/* Compute the fields radiated by the local fields. */
		for (j = 0; j < fmaconf.numbases; ++j) {
			bscenter (fmaconf.bslist[j], cen);
			for (i = 0; i < obs->count; ++i) {
				val = fsgreen (fmaconf.k0, cen, obs->locations + 3 * i);
				buf[i] += fact * val * currents[j];
			}
		}
	} else {
		/* The Green's function has already been precomputed. */
		for (j = 0, gptr = grf; j < fmaconf.numbases; ++j) {
			for (i = 0; i < obs->count; ++i, ++gptr) {
				buf[i] += fact * currents[j] * (*gptr);
			}
		}
	}


	/* Collect the solution across all processors. */
	MPI_Allreduce (buf, result, 2 * obs->count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	free (buf);

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

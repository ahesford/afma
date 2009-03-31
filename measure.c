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

int farfield (complex float *currents, measdesc *obs, complex float *result) {
	int i, mpirank;
	float *thetas, dtheta;
	Complex *fields = NULL, *crts;
	complex float fact;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* The theta samples. */
	thetas = malloc (obs->ntheta * sizeof(float));

	/* Conversion between intrinsic types and ScaleME types. */
	crts = malloc (fmaconf.numbases * sizeof(Complex)); 
	if (!mpirank) fields = malloc (obs->count * sizeof(Complex));

	dtheta = (obs->trange[1] - obs->trange[0]);
	dtheta /= MAX(obs->ntheta - 1, 1);

	for (i = 0; i < obs->ntheta; ++i)
		thetas[i] = obs->trange[0] + i * dtheta;

	for (i = 0; i < fmaconf.numbases; ++i) {
		crts[i].re = creal(currents[i]);
		crts[i].im = cimag(currents[i]);
	}

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, obs->ntheta, obs->nphi, thetas,
				obs->prange, crts, &fields))
		return 0;

	if (!mpirank) {
		/* The far-field pattern already has a factor of k in the
		 * front. However, the actual integral needs (k^2 / 4 pi), so
		 * we need the extra factors in the field. Also scale the
		 * pattern to produce actual measurements on the desired
		 * sphere. */
		fact = fmaconf.k0 * cexp (I * fmaconf.k0 * obs->radius) / (4 * M_PI * obs->radius);
		for (i = 0; i < obs->count; ++i)
			result[i] = fact * (fields[i].re + I * fields[i].im);
	}

	/* Distribute the solution to all processors. */
	MPI_Bcast (result, 2 * obs->count, MPI_FLOAT, 0, MPI_COMM_WORLD);

	free (thetas);
	free (crts);
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
			/* Don't use the radius with plane-waves. */
			rst = sin (theta);
			desc->locations[3 * k] = rst * cos (phi);
			desc->locations[3 * k + 1] = rst * sin (phi);
			desc->locations[3 * k + 2] = cos (theta);
			phi += dphi;
		}
		theta += dtheta;
	}

	return desc->count;
}

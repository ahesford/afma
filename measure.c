#include <string.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include <mpi.h>

#include "ScaleME.h"

#include "measure.h"
#include "fsgreen.h"
#include "mlfma.h"
#include "util.h"

/* Slow computation of incident field for a single source. */
int buildrhs (complex float *rhs, float *srcloc) {
	int i, j, idx[3];
	float ctr[3], off[3];
	complex float *rptr;

	for (rptr = rhs, i = 0; i < fmaconf.numbases; ++i) {
		/* Find the center of the group. */
		bscenter (fmaconf.bslist[i], off);

		/* The offset of the first basis function in the group. */
		off[0] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[1] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[2] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;

		for (j = 0; j < fmaconf.bspboxvol; ++j, ++rptr) {
			/* The position in the local grid of the basis function. */
			idx[0] = j / (fmaconf.bspbox * fmaconf.bspbox);
			idx[1] = (j / fmaconf.bspbox) % fmaconf.bspbox;
			idx[2] = j % fmaconf.bspbox;

			/* The center of the basis function. */
			ctr[0] = off[0] + fmaconf.cell * (float)idx[0];
			ctr[1] = off[1] + fmaconf.cell * (float)idx[1];
			ctr[2] = off[2] + fmaconf.cell * (float)idx[2];

			*rptr = fsplane (fmaconf.k0, ctr, srcloc) / (4 * M_PI);
		}
	}

	return 0;
}

/* Fast FMM computation of the far-field pattern using interpolation. */
int farfield (complex float *currents, measdesc *obs, complex float *result) {
	int i;
	complex float fact;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld (obs->imat, currents, &result)) return 0;

	/* The far-field pattern already has a factor of k in the front.
	 * However, the actual integral needs (k^2 / 4 pi), so we need the
	 * extra factors in the field. */
	fact = fmaconf.k0 / (4 * M_PI);
	for (i = 0; i < obs->count; ++i) result[i] *= fact;

	return 0;
}

int buildlocs (measdesc *desc) {
	int i, j, k;
	float theta, dtheta, dphi, phi, rst;

	/* Clear the interpolation matrix pointer. */
	desc->imat = NULL;

	desc->count = desc->ntheta * desc->nphi;

	dtheta = (desc->trange[1] - desc->trange[0]);
	dtheta /= MAX (desc->ntheta + 1, 1);

	dphi = (desc->prange[1] - desc->prange[0]);
	dphi /= MAX (desc->nphi, 1);

	desc->locations = malloc (3 * desc->count * sizeof(float));

	for (i = 0, k = 0; i < desc->ntheta; ++i) {
		theta = desc->trange[0] + (i + 1) * dtheta;
		for (j = 0; j < desc->nphi; ++j, ++k) {
			phi = desc->prange[0] + j * dphi;
			rst = desc->radius * sin (theta);
			desc->locations[3 * k] = rst * cos (phi);
			desc->locations[3 * k + 1] = rst * sin (phi);
			desc->locations[3 * k + 2] = desc->radius * cos (theta);
		}
	}

	return desc->count;
}

void delmeas (measdesc *desc) {
	free (desc->locations);

	/* Clear the root interpolation matrix, if applicable. */
	if (desc->imat) ScaleME_delRootInterpMat (&(desc->imat));
}

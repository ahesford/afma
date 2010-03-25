#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "ScaleME.h"

#include "mlfma.h"
#include "integrate.h"
#include "fsgreen.h"
#include "excite.h"
#include "measure.h"

#include "util.h"

/* Computes the entire RHS vector. */
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
			*rptr *= fmaconf.cellvol;
		}
	}

	return 0;
}

int multirhs (complex float *rhs, measdesc *src, complex float *mag) {
	int i, j, k, idx[3];
	float ctr[3], off[3];
	complex float *rptr;

	for (rptr = rhs, i = 0; i < fmaconf.numbases; ++i) {
		/* Find the center of the group. */
		bscenter (fmaconf.bslist[i], off);

		/* The offset of the first basis function in the group. */
		off[0] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[1] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[2] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;

		for (k = 0; k < fmaconf.bspboxvol; ++k, ++rptr) {
			/* The position in the local grid of the basis function. */
			idx[0] = k / (fmaconf.bspbox * fmaconf.bspbox);
			idx[1] = (k / fmaconf.bspbox) % fmaconf.bspbox;
			idx[2] = k % fmaconf.bspbox;

			/* The center of the basis function. */
			ctr[0] = off[0] + fmaconf.cell * (float)idx[0];
			ctr[1] = off[1] + fmaconf.cell * (float)idx[1];
			ctr[2] = off[2] + fmaconf.cell * (float)idx[2];

			*rptr = 0;

			for (j = 0; j < src->count; ++j)
				*rptr += mag[j] * fsplane (fmaconf.k0, ctr, src->locations + 3 * j) / (4 * M_PI);
			*rptr *= fmaconf.cellvol;
		}
	}

	return 0;
}

int precompgrf (measdesc *src, complex float *grf) {
	int i, j, k, idx[3];
	complex float *ptr;
	float ctr[3], off[3];

	for (ptr = grf, i = 0; i < fmaconf.numbases; ++i) {
		/* Find the center of the group. */
		bscenter (fmaconf.bslist[i], off);

		/* The offset of the first basis function in the group. */
		off[0] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[1] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[2] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;

		for (k = 0; k < fmaconf.bspboxvol; ++k) {
			/* The position in the local grid of the basis function. */
			idx[0] = k / (fmaconf.bspbox * fmaconf.bspbox);
			idx[1] = (k / fmaconf.bspbox) % fmaconf.bspbox;
			idx[2] = k % fmaconf.bspbox;

			/* The center of the basis function. */
			ctr[0] = off[0] + fmaconf.cell * (float)idx[0];
			ctr[1] = off[1] + fmaconf.cell * (float)idx[1];
			ctr[2] = off[2] + fmaconf.cell * (float)idx[2];

			for (j = 0; j < src->count; ++j, ++ptr) {
				*ptr = fsplane (fmaconf.k0, ctr, src->locations + 3 * j) / (4 * M_PI);
				*ptr *= fmaconf.cellvol;
			}
		}
	}

	return 0;
}

int precomprhs (complex float *rhs, measdesc *src, complex float *mag, complex float *grf) {
	int i, j, k;
	complex float *ptr, *rptr;

	for (i = 0, ptr = grf, rptr = rhs; i < fmaconf.numbases; ++i) {
		for (k = 0; k < fmaconf.bspboxvol; ++k, ++rptr) {
			*rptr = 0;
			for (j = 0; j < src->count; ++j, ++ptr)
				*rptr += mag[j] * (*ptr);
		}
	}

	return 0;
}

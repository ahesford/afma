#include <string.h>

#include <mpi.h>

#include "ScaleME.h"

#include "precision.h"

#include "measure.h"
#include "fsgreen.h"
#include "mlfma.h"
#include "util.h"

/* Slow computation of incident field for a single source. */
int buildrhs (cplx *rhs, real *srcloc, int plane) {
	real scale = 1.0;
	cplx (*rhsfunc) (real, real *, real *) = fsgreen;

	/* Use a plane wave instead of a point source. */
	if (plane) {
		rhsfunc = fsplane;
		scale = 1.0 / (4.0 * M_PI);
	}

#pragma omp parallel default(shared)
{
	cplx *rptr;
	real ctr[3], off[3];
	int i, j, idx[3];

#pragma omp for
	for (i = 0; i < fmaconf.numbases; ++i) {
		rptr = rhs + i * fmaconf.bspboxvol;

		/* Find the center of the group. */
		bscenter (fmaconf.bslist[i], off);

		/* The offset of the first basis function in the group. */
		off[0] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[1] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;
		off[2] += 0.5 * fmaconf.cell - 0.5 * fmaconf.grplen;

		for (j = 0; j < fmaconf.bspboxvol; ++j, ++rptr) {
			/* The position in the local grid of the basis function. */
			GRID(idx, j, fmaconf.bspbox, fmaconf.bspbox);

			/* The center of the basis function. */
			ctr[0] = off[0] + fmaconf.cell * (real)idx[0];
			ctr[1] = off[1] + fmaconf.cell * (real)idx[1];
			ctr[2] = off[2] + fmaconf.cell * (real)idx[2];

			*rptr = scale * rhsfunc (fmaconf.k0, ctr, srcloc);
		}
	}
}

	return 0;
}

/* Fast FMM computation of the far-field pattern using interpolation. */
int farfield (cplx *currents, measdesc *obs, cplx *result) {
	int i;
	cplx fact;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld (obs->imat[0], currents, &result)) return 0;

	/* The far-field pattern already has a factor of k in the front.
	 * However, the actual integral needs (k^2 / 4 pi), so we need the
	 * extra factors in the field. */
	fact = fmaconf.k0 / (4 * M_PI);
	for (i = 0; i < obs->count; ++i) result[i] *= fact;

	return 0;
}

/* Split a string str, separated by delimitors dlm, into an array of strings
 * with maximum length len. */
static int splitstr (char **arr, char *str, char *dlm, int len) {
	char **ap;

	/* Parse the string. */
	for (ap = arr; ap < arr + len; ++ap) *ap = strsep(&str, dlm);

	return len;
}

/* Attempt to an angular range for source specification. The value maxdef
 * specifies the default maximum value for the angle. */
static int measrange (char *str, real *min, real *max, real maxdef) {
	char **ap, *rn[3];
	int nr = 0;

	/* Set default angular limits values first. */
	*min = 0.0;
	*max = maxdef;

	/* Split the measurement range into an array. */
	splitstr (rn, str, ":", 3);

	/* Grab the number of angles in the range. */
	if (rn[0] && strlen(rn[0])) nr = strtol (rn[0], NULL, 0);
	/* Grab the minimum angle, if specified. */
	if (rn[1] && strlen(rn[1])) *min = strtod (rn[1], NULL);
	/* Grab the maximum angle, if specified. */
	if (rn[2] && strlen(rn[2])) *max = strtod (rn[2], NULL);

	return nr;
}

/* Build the location descriptor. */
static int buildlocs (measdesc *desc, int plane, int ntheta,
		real tmin, real tmax, int nphi, real pmin, real pmax) {
	int i, j, k;
	real theta, dtheta, dphi, phi, rst;

	/* Clear the interpolation matrix pointer. */
	desc->imat[0] = desc->imat[1] = NULL;

	/* Configure the locations to be plane-wave directions or point sources. */
	desc->plane = plane;

	/* Record the number of angular samples desired. */
	desc->ntheta = ntheta;
	desc->nphi = nphi;

	/* Record the limits of the ranges. Convert to radians. */
	desc->trange[0] = tmin * M_PI / 180.0;
	desc->trange[1] = tmax * M_PI / 180.0;
	desc->prange[0] = pmin * M_PI / 180.0;
	desc->prange[1] = pmax * M_PI / 180.0;

	/* Count the total number of measurements and allocate the location array. */
	desc->count = desc->ntheta * desc->nphi;
	desc->locations = malloc (3 * desc->count * sizeof(real));

	/* Calculate the angular steps. The polar angle avoids the poles. */
	dtheta = (desc->trange[1] - desc->trange[0]);
	dtheta /= MAX (desc->ntheta + 1, 1);
	dphi = (desc->prange[1] - desc->prange[0]);
	dphi /= MAX (desc->nphi, 1);

	/* Populate the location array. Skip the north pole! */
	for (i = 0, k = 0; i < desc->ntheta; ++i) {
		theta = desc->trange[0] + (i + 1) * dtheta;
		for (j = 0; j < desc->nphi; ++j, ++k) {
			phi = desc->prange[0] + j * dphi;
			rst = sin (theta);
			desc->locations[3 * k] = rst * cos (phi);
			desc->locations[3 * k + 1] = rst * sin (phi);
			desc->locations[3 * k + 2] = cos (theta);
		}
	}

	return desc->count;
}

/* Build the source descriptor with specification string spec. */
int buildsrc (measdesc *desc, char *spec) {
	real tmin = 0, tmax = 0, pmin = 0, pmax = 0, srcpt[3];
	char *srcloc[3];
	int ntheta = 0, nphi = 0, plane = 1;

	/* Split the source locations into a theta and phi range, or coordinates. */
	splitstr (srcloc, spec, ",", 3);

	if (srcloc[2]) {
		/* If a third argument is present, use a point source. */
		plane = 0;
		srcpt[0] = strtod (srcloc[0], NULL);
		srcpt[1] = strtod (srcloc[1], NULL);
		srcpt[2] = strtod (srcloc[2], NULL);

		ntheta = nphi = 1;
	} else {
		/* Otherwise, use a range of plane waves. */
		ntheta = measrange (srcloc[0], &tmin, &tmax, 180);
		nphi = measrange (srcloc[1], &pmin, &pmax, 360);
	}

	/* Build the observer locations. */
	buildlocs (desc, plane, ntheta, tmin, tmax, nphi, pmin, pmax);

	/* Override the location for a point source. */
	if (!plane) memcpy (desc->locations, srcpt, 3 * sizeof(real));

	return desc->count;
}

/* Build the observation descriptor with specification string spec. */
int buildobs (measdesc *desc, char *spec) {
	real tmin = 0, tmax = 0, pmin = 0, pmax = 0;
	char *obsloc[2];
	int ntheta = 0, nphi = 0;

	/* Split the observer locations into a theta and phi range. */
	splitstr (obsloc, spec, ",", 2);

	/* Parse the theta and phi ranges. */
	ntheta = measrange (obsloc[0], &tmin, &tmax, 180);
	nphi = measrange (obsloc[1], &pmin, &pmax, 360);

	/* Build the observer locations. The type is always a plane-wave. */
	buildlocs (desc, 1, ntheta, tmin, tmax, nphi, pmin, pmax);

	return desc->count;
}

void delmeas (measdesc *desc) {
	free (desc->locations);

	/* Clear the root interpolation matrix, if applicable. */
	if (desc->imat[0]) ScaleME_delRootInterpMat (desc->imat);
	if (desc->imat[1]) ScaleME_delRootInterpMat (desc->imat + 1);
}

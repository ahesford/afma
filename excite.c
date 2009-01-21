#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "mlfma.h"
#include "integrate.h"
#include "fsgreen.h"
#include "excite.h"
#include "measure.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Computes the RHS for a plane-wave in a given direction. */
complex float planerhs (int gi, float *srcdir) {
	complex float ans;
	float ctr[3];

	/* Find the center of the requested basis. */
	bscenter (gi, ctr);

	/* Use a single-point integration rule. */
	ans = fsplane (fmaconf.k0, ctr, srcdir);
	ans *= fmaconf.cellvol;

	return ans;
}

/* Computes the RHS for a point source at the given location. */
complex float pointrhs (int gi, float *srcloc) {
	complex float ans;
	float ctr[3];

	bscenter (gi, ctr);
	ans = rcvint (fsgreen, fmaconf.k0, srcloc, ctr, fmaconf.cell,
			fmaconf.nrcvint, fmaconf.rcvpts, fmaconf.rcvwts);

	return ans;
}

/* Computes the entire RHS vector. */
int buildrhs (complex float *rhs, float *srcloc, int type) {
	int i;
	complex float (*rhsfunc) (int, float *);

	if (type) rhsfunc = planerhs;
	else rhsfunc = pointrhs;

	for (i = 0; i < fmaconf.numbases; ++i)
		rhs[i] = rhsfunc (fmaconf.bslist[i], srcloc);

	return fmaconf.numbases;
}

int multirhs (complex float *rhs, measdesc *src, complex float *mag, int type) {
	int i, j;
	complex float (*rhsfunc) (int, float *);

	if (type) rhsfunc = planerhs;
	else rhsfunc = pointrhs;

	for (i = 0; i < fmaconf.numbases; ++i) {
		rhs[i] = 0;
		for (j = 0; j < src->count; ++j)
			rhs[i] += mag[j] * rhsfunc (fmaconf.bslist[i], src->locations + 3 * j);
	}

	return fmaconf.numbases * src->count;
}

int precompgrf (measdesc *src, complex float *grf, int type) {
	int i, j;
	complex float *ptr = grf;
	complex float (*rhsfunc) (int, float *);

	if (type) rhsfunc = planerhs;
	else rhsfunc = pointrhs;

	for (i = 0; i < fmaconf.numbases; ++i) {
		for (j = 0; j < src->count; ++j) {
			*ptr = rhsfunc (fmaconf.bslist[i], src->locations + 3 * j);
			++ptr;
		}
	}

	return fmaconf.numbases * src->count;
}

int precomprhs (complex float *rhs, measdesc *src, complex float *mag, complex float *grf) {
	int i, j;
	complex float *ptr = grf;

	for (i = 0; i < fmaconf.numbases; ++i) {
		rhs[i] = 0;
		for (j = 0; j < src->count; ++j) {
			rhs[i] += mag[j] * (*ptr);
			++ptr;
		}
	}

	return fmaconf.numbases * src->count;
}

int farexcite (complex float *excit, measdesc *obs, complex float *currents) {
	int i;
	float *thetas, dtheta;
	Complex *fields = (Complex *)excit;

	thetas = malloc (obs->ntheta * sizeof(float));

	dtheta = (obs->trange[1] - obs->trange[0]);
	dtheta /= MAX(obs->ntheta - 1, 1);

	for (i = 0; i < obs->ntheta; ++i)
		thetas[i] = obs->trange[0] + i * dtheta;

	/* The result must be a pointer to the pointer. */
	if (ScaleME_setRootFarFld_PI (6, obs->ntheta, obs->nphi, thetas,
				obs->prange, (Complex *)currents, &fields))
		return 0;

	free (thetas);
	return 1;
}

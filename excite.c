#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "mlfma.h"
#include "integrate.h"
#include "fsgreen.h"
#include "excite.h"
#include "measure.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Computes the entire RHS vector. */
int buildrhs (complex float *rhs, float *srcloc, int type) {
	complex float (*grf) (float, float *, float *);
	float ctr[3];
	int i;

	if (type) grf = fsplane;
	else grf = fsgreen;

	for (i = 0; i < fmaconf.numbases; ++i) {
		bscenter (fmaconf.bslist[i], ctr);
		rhs[i] = grf (fmaconf.k0, ctr, srcloc);
	}

	return fmaconf.numbases;
}

int multirhs (complex float *rhs, measdesc *src, complex float *mag, int type) {
	int i, j;
	float ctr[3];
	complex float (*grf) (float, float *, float *);

	if (type) grf = fsplane;
	else grf = fsgreen;

	for (i = 0; i < fmaconf.numbases; ++i) {
		rhs[i] = 0;
		bscenter (fmaconf.bslist[i], ctr);
		for (j = 0; j < src->count; ++j)
			rhs[i] += mag[j] * grf (fmaconf.k0, ctr, src->locations + 3 * j);
	}

	return fmaconf.numbases * src->count;
}

int precompgrf (measdesc *src, complex float *grf, int type) {
	int i, j;
	float ctr[3];
	complex float *ptr = grf;
	complex float (*rhsfunc) (float, float *, float *);

	if (type) rhsfunc = fsplane;
	else rhsfunc = fsgreen;

	for (i = 0; i < fmaconf.numbases; ++i) {
		bscenter (fmaconf.bslist[i], ctr);
		for (j = 0; j < src->count; ++j) {
			*ptr = rhsfunc (fmaconf.k0, ctr, src->locations + 3 * j);
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

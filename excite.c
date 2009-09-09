#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "ScaleME.h"

#include "mlfma.h"
#include "integrate.h"
#include "fsgreen.h"
#include "excite.h"
#include "measure.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Computes the entire RHS vector. */
int buildrhs (complex float *rhs, float *srcloc) {
	int i;
	float ctr[3];

	for (i = 0; i < fmaconf.numbases; ++i) {
		bscenter (fmaconf.bslist[i], ctr);
		rhs[i] = fsplane (fmaconf.k0, ctr, srcloc);
	}

	return fmaconf.numbases;
}

int multirhs (complex float *rhs, measdesc *src, complex float *mag) {
	int i, j;
	float ctr[3];

	for (i = 0; i < fmaconf.numbases; ++i) {
		rhs[i] = 0;
		bscenter (fmaconf.bslist[i], ctr);
		for (j = 0; j < src->count; ++j)
			rhs[i] += mag[j] * fsplane (fmaconf.k0, ctr, src->locations + 3 * j);
	}

	return fmaconf.numbases * src->count;
}

int precompgrf (measdesc *src, complex float *grf) {
	int i, j;
	complex float *ptr = grf;
	float ctr[3];

	for (i = 0; i < fmaconf.numbases; ++i) {
		bscenter (fmaconf.bslist[i], ctr);
		for (j = 0; j < src->count; ++j) {
			*ptr = fsplane (fmaconf.k0, ctr, src->locations + 3 * j);
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

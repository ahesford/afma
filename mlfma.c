#include <complex.h> // This is provided by GCC.
#include <Complex.h> // This is provided by ScaleME.

#include "utility.h"
#include "integrate.h"
#include "mlfma.h"

fmadesc mlfma;

/* Compute the radiation pattern of a square cell. */
void radpattern (int gi, float *cen, float *s, Complex *ans) {
	float sc, sr, rv[3];
	complex float val;

	bscenter (gi, rv);

	sc = s[0] * cen[0] + s[1] * cen[1] + s[2] * cen[2];
	sr = s[0] * rv[0] + s[1] * rv[1] + s[2] * rv[2];

	val = mlfma.k0 * cexp (I * mlfma.k0 * sc) * cexp (-I * mlfma.k0 * sr);
	val *= mlfma.cell[0] * mlfma.cell[1] * mlfma.cell[2];

	/* Copy the value into the solution. */
	ans->re = creal (val);
	ans->im = cimag (val);
}

/* The receiving pattern is just the conjugate of the radiation pattern. */
void rcvpattern (int gi, float *cen, float *s, Complex *ans) {
	radpattern (gi, cen, s, ans);
	ans->im = -(ans->im);
}

/* Direct interactions between two global basis indices. */
void impedance (int gi, int gj, Complex *ans) {
	float src[3], obs[3];
	complex float val;

	/* This is the self term, so use the analytic approximation. */
	if (gi == gj) val = selfint (mlfma.k0, mlfma.cell);
	else { 
		/* Find the source and observation centers. */
		bscenter (gi, obs);
		bscenter (gj, src); 
		
		/* Compute the near-neighbor term, using Gaussian quadrature. */
		val = srcint (mlfma.k0, src, obs, mlfma.cell);
	}
	
	/* Copy the value into the solution. */
	ans->re = creal (val);
	ans->im = cimag (val);
	return;
}

/* Find the center of the basis function referred to by the global index. */
void bscenter (int gi, float *cen) {
	int i, j, k;

	i = gi % mlfma.nx;
	j = (gi / mlfma.nx) % mlfma.ny;
	k = gi / (mlfma.nx * mlfma.ny);

	cen[0] = mlfma.min[0] + ((float)i + 0.5) * mlfma.cell[0];
	cen[1] = mlfma.min[1] + ((float)j + 0.5) * mlfma.cell[1];
	cen[2] = mlfma.min[2] + ((float)k + 0.5) * mlfma.cell[2];
}

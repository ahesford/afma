#include <complex.h> // This is provided by GCC.
#include <math.h>

#include <Complex.h> // This is provided by ScaleME.
#include <ScaleME.h>

#include "fsgreen.h"
#include "integrate.h"
#include "mlfma.h"

fmadesc fmaconf;

static complex float genpat (int gi, float *cen, float *s) {
	float sc, sr, rv[3];
	complex float val;

	bscenter (gi, rv);

	sc = s[0] * cen[0] + s[1] * cen[1] + s[2] * cen[2];
	sr = s[0] * rv[0] + s[1] * rv[1] + s[2] * rv[2];

	val = cexp (I * fmaconf.k0 * sc) * cexp (-I * fmaconf.k0 * sr);
	val *= fmaconf.cellvol;

	return val;
}

/* Compute the radiation pattern of a square cell. The formulation assumes
 * the Green's function does not have a 4 pi term in the denominator. */
void radpattern (int gi, float *cen, float *s, Complex *ans) {
	complex float val;

	/* The k0 term is due to the limit of the radiation pattern. */
	val = fmaconf.k0 * genpat (gi, cen, s);

	/* Copy the value into the solution. */
	ans->re = creal (val);
	ans->im = cimag (val);
}

/* The receiving pattern is the conjugate of the radiation pattern. */
void rcvpattern (int gi, float *cen, float *s, Complex *ans) {
	complex float val;

	val = conj (genpat (gi, cen, s));
	/* Include constants outside the integral. The 4 pi factor is included
	 * because it was neglected in radiation pattern construction. */
	val *= fmaconf.k0 * fmaconf.k0 / (4 * M_PI);

	ans->re = creal (val);
	ans->im = cimag (val);
}

/* Direct interactions between two global basis indices. */
void impedance (int gi, int gj, Complex *ans) {
	float obs[3], src[3];
	complex float val;

	/* This is the self term, so use the analytic approximation. */
	if (gi == gj) val = selfint (fmaconf.k0, fmaconf.cell);
	else {
		/* Find the source and observation centers. */
		bscenter (gi, obs);
		bscenter (gj, src);

		val = srcint (fmaconf.k0, src, obs, fmaconf.cell,
				fmaconf.nsrcint, fmaconf.srcpts, fmaconf.srcwts,
				fmaconf.nrcvint, fmaconf.rcvpts, fmaconf.rcvwts);
		val *= fmaconf.k0 * fmaconf.k0;
	}

	/* Copy the value into the solution. */
	ans->re = creal (val);
	ans->im = cimag (val);
	return;
}

void bscenter (int gi, float *cen) {
	int i, j, k;

	i = gi % fmaconf.nx;
	j = (gi / fmaconf.nx) % fmaconf.ny;
	k = gi / (fmaconf.nx * fmaconf.ny);

	cen[0] = fmaconf.min[0] + ((float)i + 0.5) * fmaconf.cell[0];
	cen[1] = fmaconf.min[1] + ((float)j + 0.5) * fmaconf.cell[1];
	cen[2] = fmaconf.min[2] + ((float)k + 0.5) * fmaconf.cell[2];
}

void bsindex (int gi, int *idx) {
	idx[0] = gi % fmaconf.nx;
	idx[1] = (gi / fmaconf.nx) % fmaconf.ny;
	idx[2] = gi / (fmaconf.nx * fmaconf.ny);
}

void interaction (int gi, int gj, Complex *ans) {
	int idx, dist[3], obs[3], src[3];

	bsindex (gj, src);
	bsindex (gi, obs);

	dist[0] = abs (src[0] - obs[0]);
	dist[1] = abs (src[1] - obs[1]);
	dist[2] = abs (src[2] - obs[2]);

	if (dist[0] >= fmaconf.nbors[0] || dist[1] >= fmaconf.nbors[1]
			|| dist[2] >= fmaconf.nbors[2]) {
		fprintf (stderr, "Interaction (%d -> %d) not precomputed\n", gj, gi);
		impedance (gi, gj, ans);
		return;
	}

	/* Find the linear index for the interaction. */
	idx = dist[2] + dist[1] * fmaconf.nbors[2]
		+ dist[0] * fmaconf.nbors[2] * fmaconf.nbors[1];

	/* Copy the value into the solution. */
	ans->re = creal (fmaconf.gridints[idx]);
	ans->im = cimag (fmaconf.gridints[idx]);

	return;
}

int preimpedance () {
	float lbox[3], dist[3], zero[3] = { 0, 0, 0 };
	int i, j, k, idx, totel;

	/* Find the smallest box size. */
	ScaleME_getSmallestBoxSize (lbox);

	/* Find the maximum number of near neighbors in each dimension. */
	fmaconf.nbors[0] = (int) (((fmaconf.numbuffer + 1) * lbox[0]) / fmaconf.cell[0]) + 1;
	fmaconf.nbors[1] = (int) (((fmaconf.numbuffer + 1) * lbox[1]) / fmaconf.cell[1]) + 1;
	fmaconf.nbors[2] = (int) (((fmaconf.numbuffer + 1) * lbox[2]) / fmaconf.cell[2]) + 1;

	fprintf (stderr, "Precomputing interactions for (%d, %d, %d) grid\n",
			fmaconf.nbors[0], fmaconf.nbors[1], fmaconf.nbors[2]);

	totel = fmaconf.nbors[0] * fmaconf.nbors[1] * fmaconf.nbors[2];

	fmaconf.gridints = malloc (totel * sizeof (complex float));

	/* Compute the interactions. */
	for (i = 0, idx = 0; i < fmaconf.nbors[0]; ++i) {
		dist[0] = i * fmaconf.cell[0];
		for (j = 0; j < fmaconf.nbors[1]; ++j) {
			dist[1] = j * fmaconf.cell[1];
			for (k = 0; k < fmaconf.nbors[2]; ++k, ++idx) {
				dist[2] = k * fmaconf.cell[2];
				
				fmaconf.gridints[idx] = srcint (fmaconf.k0, zero, dist, fmaconf.cell,
						fmaconf.nsrcint, fmaconf.srcpts, fmaconf.srcwts,
						fmaconf.nrcvint, fmaconf.rcvpts, fmaconf.rcvwts);
				fmaconf.gridints[idx] *= fmaconf.k0 * fmaconf.k0;
			}
		}
	}

	/* Recompute the self interaction with high accuracy. */
	fmaconf.gridints[0] = selfinthigh (fmaconf.k0, fmaconf.cell);

	return totel;
}

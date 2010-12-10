#include "precision.h"

#include "fsgreen.h"

/* Computes the free-space Green's function between two points. */
cplx fsgreen (real k, real *r, real *rp) {
	cplx ans;
	real dist, kdist;

	/* Compute the distance between elements. */
	dist = (r[0] - rp[0]) * (r[0] - rp[0]);
	dist += (r[1] - rp[1]) * (r[1] - rp[1]);
	dist += (r[2] - rp[2]) * (r[2] - rp[2]);
	dist = sqrt (dist);

	kdist = k * dist;

	ans = (cos (kdist) + I * sin (kdist))  / (4 * M_PI * dist);

	return ans;
}

/* Computes the free-space Green's function between two points,
 * with the singular part subtracted off to improve integration accuracy. */
cplx fsgrnsmooth (real k, real *r, real *rp) {
	real dist;

	/* Compute the distance between elements. */
	dist = (r[0] - rp[0]) * (r[0] - rp[0]);
	dist += (r[1] - rp[1]) * (r[1] - rp[1]);
	dist += (r[2] - rp[2]) * (r[2] - rp[2]);
	dist = sqrt (dist);

	/* For zero values, evaluate the limit analytically. */
	if (dist < REAL_EPSILON) return I * k / (4 * M_PI);

	/* Compute the smooth part numerically. */
	return fsgreen (k, r, rp) - 1. / (4 * M_PI * dist);
}

/* Computes a plane wave from a specific direction at a point. */
cplx fsplane (real k, real *r, real *s) {
	real sr, ds, ksr;

	ds = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
	sr = (s[0] * r[0] + s[1] * r[1] + s[2] * r[2]) / ds;

	ksr = k * sr;

	return cos (ksr) - I * sin (ksr);
}

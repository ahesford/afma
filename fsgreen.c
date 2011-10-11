#include "precision.h"

#include "fsgreen.h"

/* Computes the free-space Green's function between two points. */
cplx fsgreen (real k, real *r, real *rp) {
	cplx ans;
	real dist;

	/* Compute the distance between elements. */
	dist = (r[0] - rp[0]) * (r[0] - rp[0]);
	dist += (r[1] - rp[1]) * (r[1] - rp[1]);
	dist += (r[2] - rp[2]) * (r[2] - rp[2]);
	dist = sqrt (dist);

	ans = cexp(I * k * dist) / (4 * M_PI * dist);

	return ans;
}

/* Computes the free-space Green's function with the source in coordinates
 * transformed according to Duffy's rule and the observation at the origin. The
 * argument rv is ignored but is present so the function can be plugged into
 * the source integration routine. */
cplx fsgrnduffy (real k, real *r, real *rv) {
	cplx ans;
	real dist;

	/* Compute the distance between elements. */
	dist = sqrt(1. + r[1] * r[1] + r[2] * r[2]);

	ans = r[0] * cexp(I * k * r[0] * dist) / (4 * M_PI * dist);

	return ans;
}

/* Computes a plane wave from a specific direction at a point. */
cplx fsplane (real k, real *r, real *s) {
	real sr, ds;

	ds = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
	sr = (s[0] * r[0] + s[1] * r[1] + s[2] * r[2]) / ds;

	return cexp(-I * k * sr);
}

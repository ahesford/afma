#include <math.h>
#include <complex.h>

#include "fsgreen.h"

/* Computes the free-space Green's function between two points. */
complex float fsgreen (float k, float *r, float *rp) {
	complex float ans;
	float dist;

	/* Compute the distance between elements. */
	dist = (r[0] - rp[0]) * (r[0] - rp[0]);
	dist += (r[1] - rp[1]) * (r[1] - rp[1]);
	dist += (r[2] - rp[2]) * (r[2] - rp[2]);
	dist = sqrt (dist);

	ans = cexp (I * k * dist) / (4 * M_PI * dist);

	return ans;
}

/* Computes a plane wave from a specific direction at a point. */
complex float fsplane (float k, float *r, float *s) {
	complex float ans;
	float sr;

	sr = s[0] * r[0] + s[1] * r[1] + s[2] * r[2];

	ans = cexp (-I * k * sr) / (4 * M_PI);

	return ans;
}

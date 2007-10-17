#include <math.h>
#include <complex.h>

complex float fsgreen (complex float k, float *r, float *rp) {
	complex float ans;
	float dist;

	/* Compute the distance between elements. */
	dist = (r[0] - rp[0]) * (r[0] - rp[0]);
	dist += (r[1] - rp[1]) * (r[1] - rp[1]);
	dist += (r[2] - rp[2]) * (r[2] - rp[2]);
	dist = sqrt (dist);

	ans = cexp (I * k * dist) / dist;

	return ans;
}

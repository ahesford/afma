#include <math.h>
#include <complex.h>

#include "fsgreen.h"
#include "integrate.h"

float pts[4] = { -OUTPT, -INPT, INPT, OUTPT };
float wts[4] = { OUTWT, INWT, INWT, OUTWT };

/* Four-point (per dimension) integration of the source. */
complex float srcint (float k, float *src, float *obs, float *dc) {
	complex float ans = 0, val;
	int i, j, l;
	float srcpt[3];

	for (i = 0; i < NUMPTS; ++i) {
		srcpt[0] = src[0] + 0.5 * dc[0] * pts[i];
		for (j = 0; j < NUMPTS; ++j) {
			srcpt[1] = src[1] + 0.5 * dc[1] * pts[j];
			for (l = 0; l < NUMPTS; ++l) {
				srcpt[2] = src[2] + 0.5 * dc[2] * pts[l];
				val = fsgreen (k, srcpt, obs);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc[0] * dc[1] * dc[2] / 8;
	return ans;
}

/* Analytic approximation to the self-integration term. The square cell is
 * approximated as a sphere with the same volume. */
complex float selfint (float k, float *dc) {
	complex float ans, ikr;
	float r;

	r = cbrt (3 * dc[0] * dc[1] * dc[2] / (4 * M_PI));
	ikr = I * k * r;

	ans = (1.0 - ikr) * cexp (ikr) - 1.0;

	return ans;
}

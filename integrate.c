#include <math.h>
#include <complex.h>

#include "utility.h"
#include "integrate.h"

float pts[4] = { -OUTPT, -INPT, INPT, OUTPT };
float wts[4] = { OUTWT, INWT, INWT, OUTWT };

complex float rcvint (igrandf green, float k, float *src, float *cen, float *dc) {
	complex float ans = 0, val;
	int i, j, l;
	float obs[3];

	for (i = 0; i < NUMPTS; ++i) {
		obs[0] = cen[0] + 0.5 * dc[0] * pts[i];
		for (j = 0; j < NUMPTS; ++j) {
			obs[1] = cen[1] + 0.5 * dc[1] * pts[j];
			for (l = 0; l < NUMPTS; ++l) {
				obs[2] = cen[2] + 0.5 * dc[2] * pts[l];
				val = green (k, obs, src);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc[0] * dc[1] * dc[2] / 8;
	return ans;
}

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
				val = rcvint (fsgreen, k, srcpt, obs, dc);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc[0] * dc[1] * dc[2] / 8;
	return ans;
}

complex float selfint (float k, float *dc) {
	complex float ans, ikr;
	float r;

	r = cbrt (3 * dc[0] * dc[1] * dc[2] / (4 * M_PI));
	ikr = I * k * r;

	ans = (1.0 - ikr) * cexp (ikr) - 1.0;
	ans *= (4 * M_PI) / (k * k);

	return ans;
}

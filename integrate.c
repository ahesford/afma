#include <math.h>
#include <complex.h>

#include "fsgreen.h"
#include "integrate.h"

/* Three-point (per dimension) integration of the observer. */
complex float rcvint (float k, float *src, float *obs, float dc, ifunc grf) {
	const int numpts = 3;
	complex float ans = 0, val;
	int i, j, l;
	float obspt[3];
	float pts[3] = { -OTPT3, 0.00, OTPT3 };
	float wts[3] = { OTWT3, CNWT3, OTWT3 };

	for (i = 0; i < numpts; ++i) {
		obspt[0] = obs[0] + 0.5 * dc * pts[i];
		for (j = 0; j < numpts; ++j) {
			obspt[1] = obs[1] + 0.5 * dc * pts[j];
			for (l = 0; l < numpts; ++l) {
				obspt[2] = obs[2] + 0.5 * dc * pts[l];
				val = srcint (k, src, obspt, dc, grf);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc * dc * dc / 8;
	return ans;
}

/* Four-point (per dimension) integration of the source. */
complex float srcint (float k, float *src, float *obs, float dc, ifunc grf) {
	const int numpts = 4;
	complex float ans = 0, val;
	int i, j, l;
	float srcpt[3];
	float pts[4] = { -OTPT4, -INPT4, INPT4, OTPT4 };
	float wts[4] = { OTWT4, INWT4, INWT4, OTWT4 };


	for (i = 0; i < numpts; ++i) {
		srcpt[0] = src[0] + 0.5 * dc * pts[i];
		for (j = 0; j < numpts; ++j) {
			srcpt[1] = src[1] + 0.5 * dc * pts[j];
			for (l = 0; l < numpts; ++l) {
				srcpt[2] = src[2] + 0.5 * dc * pts[l];
				val = grf (k, srcpt, obs);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc * dc * dc / 8;
	return ans;
}

/* Analytic approximation to the self-integration term. The square cell is
 * approximated as a sphere with the same volume. */
complex float selfint (float k, float dc) {
	complex float ans, ikr;
	float r;

	r = cbrt (3 * dc * dc * dc / (4 * M_PI));
	ikr = I * k * r;

	ans = (1.0 - ikr) * cexp (ikr) - 1.0;

	return ans;
}

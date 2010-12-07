#include <stdlib.h>

#include <math.h>
#include <complex.h>

#include "fsgreen.h"
#include "integrate.h"
#include "util.h"

static float *rcvpts = NULL, *rcvwts = NULL, *srcpts = NULL, *srcwts = NULL;
static int numrcvpts = 0, numsrcpts = 0;

/* Three-point (per dimension) integration of the observer. */
complex float rcvint (float k, float *src, float *obs, float dc, ifunc grf) {
	complex float ans = 0, val;
	int i, j, l;
	float obspt[3];

	for (i = 0; i < numrcvpts; ++i) {
		obspt[0] = obs[0] + 0.5 * dc * rcvpts[i];
		for (j = 0; j < numrcvpts; ++j) {
			obspt[1] = obs[1] + 0.5 * dc * rcvpts[j];
			for (l = 0; l < numrcvpts; ++l) {
				obspt[2] = obs[2] + 0.5 * dc * rcvpts[l];
				val = srcint (k, src, obspt, dc, grf);
				ans += rcvwts[i] * rcvwts[j] * rcvwts[l] * val;
			}
		}
	}

	ans *= dc * dc * dc / 8;
	return ans;
}

/* Four-point (per dimension) integration of the source. */
complex float srcint (float k, float *src, float *obs, float dc, ifunc grf) {
	complex float ans = 0, val;
	int i, j, l;
	float srcpt[3];

	for (i = 0; i < numsrcpts; ++i) {
		srcpt[0] = src[0] + 0.5 * dc * srcpts[i];
		for (j = 0; j < numsrcpts; ++j) {
			srcpt[1] = src[1] + 0.5 * dc * srcpts[j];
			for (l = 0; l < numsrcpts; ++l) {
				srcpt[2] = src[2] + 0.5 * dc * srcpts[l];
				val = grf (k, srcpt, obs);
				ans += srcwts[i] * srcwts[j] * srcwts[l] * val;
			}
		}
	}

	ans *= dc * dc * dc / 8;
	return ans;
}

/* Use either the singularity-extracted approximation to the self-integration
 * term, with four-point source integration; or use an analytic approximation. */
complex float selfint (float k, float dc, int analytic) {
	complex float ikr;
	float r, zero[3] = {0., 0., 0.};

	r = cbrt (3. / (4. * M_PI)) * dc;

	/* Sum the contributions of the smooth and singular parts. */
	if (!analytic) return srcint (k, zero, zero, dc, fsgrnsmooth) + 0.5 * r * r;

	/* Otherwise use the analytic approximation. */
	ikr = I * k * r;
	return ((1.0 - ikr) * cexp (ikr) - 1.0) / (k * k);
}

void bldintrules (int nspts, int nrpts) {
	if (nspts > 0) {
		srcpts = malloc (2 * nspts * sizeof(float));
		srcwts = srcpts + nspts;
		gaussleg (srcpts, srcwts, nspts);
		numsrcpts = nspts;
	}
	if (nrpts > 0) {
		rcvpts = malloc (2 * nrpts * sizeof(float));
		rcvwts = rcvpts + nrpts;
		gaussleg (rcvpts, rcvwts, nrpts);
		numrcvpts = nrpts;
	}
}

void delintrules () {
	numrcvpts = numsrcpts = 0;
	if (srcpts) free (srcpts);
	if (rcvpts) free (rcvpts);

	srcpts = srcwts = rcvpts = rcvwts = NULL;
}

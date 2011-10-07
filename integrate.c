#include <stdlib.h>

#include "precision.h"

#include "fsgreen.h"
#include "integrate.h"
#include "util.h"

static real *rcvpts = NULL, *rcvwts = NULL, *srcpts = NULL, *srcwts = NULL;
static int numrcvpts = 0, numsrcpts = 0;

/* N-point (per dimension) integration of the observer. */
cplx rcvint (real k, real *src, real *obs, real dc, ifunc grf) {
	cplx ans = 0, val;
	int i, j, l;
	real obspt[3];

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

/* N-point (per dimension) integration of the source. */
cplx srcint (real k, real *src, real *obs, real dc, ifunc grf) {
	cplx ans = 0, val;
	int i, j, l;
	real spt[3];

	for (i = 0; i < numsrcpts; ++i) {
		spt[0] = src[0] + 0.5 * dc * srcpts[i];
		for (j = 0; j < numsrcpts; ++j) {
			spt[1] = src[1] + 0.5 * dc * srcpts[j];
			for (l = 0; l < numsrcpts; ++l) {
				spt[2] = src[2] + 0.5 * dc * srcpts[l];
				val = grf (k, spt, obs);
				ans += srcwts[i] * srcwts[j] * srcwts[l] * val;
			}
		}
	}

	ans *= dc * dc * dc / 8;
	return ans;
}

/* N-point (per dimension) Duffy integration of the self term. */
cplx duffyint (real k, real dc) {
	cplx ans = 0, val;
	int i, j, l;
	real spt[3];

	for (i = 0; i < numsrcpts; ++i) {
		spt[0] = 0.25 * dc * (srcpts[i] + 1.);
		for (j = 0; j < numsrcpts; ++j) {
			spt[1] = srcpts[j];
			for (l = 0; l < numsrcpts; ++l) {
				spt[2] = srcpts[l];
				val = fsgrnduffy (k, spt);
				ans += srcwts[i] * srcwts[j] * srcwts[l] * val;
			}
		}
	}

	 /* The pyramidal integration must be scaled by dc/4 due to the change
	  * of variables. This is then scaled by 6 because a pyramid is 1/6 the
	  * total cubic cell. */
	return 1.5 * dc * ans;
}

/* Use an analytic approximation to the self-integration. */
cplx selfint (real k, real dc) {
	cplx ikr;
	real r, zero[3] = {0., 0., 0.};

	r = cbrt (3. / (4. * M_PI)) * dc;

	ikr = I * k * r;
	return ((1.0 - ikr) * cexp (ikr) - 1.0) / (k * k);
}

void bldintrules (int nspts, int nrpts) {
	if (nspts > 0) {
		srcpts = malloc (2 * nspts * sizeof(real));
		srcwts = srcpts + nspts;
		gaussleg (srcpts, srcwts, nspts);
		numsrcpts = nspts;
	}
	if (nrpts > 0) {
		rcvpts = malloc (2 * nrpts * sizeof(real));
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

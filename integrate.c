#include <stdlib.h>

#include "precision.h"

#include "fsgreen.h"
#include "integrate.h"
#include "util.h"

static real *rcvpts = NULL, *rcvwts = NULL, *srcpts = NULL, *srcwts = NULL;
static int numrcvpts = 0, numsrcpts = 0;

/* Perform a cyclic rotation of the elements of the 3-D vector x. */
static inline void rotate (real *x) {
	real tmp;
	
	tmp = x[0];
	x[0] = x[1];
	x[1] = x[2];
	x[2] = tmp;

}

/* Halve the cell dimensions for easy reference in the integration routines. */
static inline void halfcell (real *hcell, real *dc) {
	hcell[0] = 0.5 * dc[0];
	hcell[1] = 0.5 * dc[1];
	hcell[2] = 0.5 * dc[2];
}

/* Compute the Duffy-transformed center and dimensions of a pyramidal portion
 * of a cell with Cartesian dimensions d and a singularity at s. The variable
 * sign governs whether the axis of the pyramid is along the positive (+1) or
 * negative (-1) x axis. */
static inline void duffcell (real *ctr, real *dc, real *s, real *d, int sign) {
	/* Force the sign to have unity magnitude. */
	sign = (sign > 0) ? 1 : -1;
	/* Compute the square-limit length. */
	real l = 0.5 * d[0] - sign * s[0];
	/* Popuplate the Duffy cell dimensions. */
	dc[0] = l;
	dc[1] = d[1] / l;
	dc[2] = d[2] / l;
	/* Populate the Duffy cell center. */
	ctr[0] = 0.5 * l;
	ctr[1] = -s[1] / l;
	ctr[2] = -s[2] / l;
}

/* N-point (per dimension) double integration over cells with dimensions dc and
 * centered on src and obs. The integrator integ can be any single integration
 * routine that acts on the integrand grf from src to obs with wave number k. */
cplx rcvint (real k, real *src, real *obs, real *dc, integrator integ, ifunc grf) {
	cplx ans = 0, val;
	int i, j, l;
	real obspt[3], hcell[3];

	halfcell(hcell, dc);

	for (i = 0; i < numrcvpts; ++i) {
		obspt[0] = obs[0] + hcell[0] * rcvpts[i];
		for (j = 0; j < numrcvpts; ++j) {
			obspt[1] = obs[1] + hcell[1] * rcvpts[j];
			for (l = 0; l < numrcvpts; ++l) {
				obspt[2] = obs[2] + hcell[2] * rcvpts[l];
				val = integ (k, src, obspt, dc, grf);
				ans += rcvwts[i] * rcvwts[j] * rcvwts[l] * val;
			}
		}
	}

	ans *= hcell[0] * hcell[1] * hcell[2];
	return ans;
}

/* N-point (per dimension) integration over a cell with dimensions dc and
 * center src of the integrand Green's function grf operating from src to obs
 * with wave number k. */
cplx srcint (real k, real *src, real *obs, real *dc, ifunc grf) {
	cplx ans = 0, val;
	int i, j, l;
	real spt[3], hcell[3];

	halfcell(hcell, dc);

	for (i = 0; i < numsrcpts; ++i) {
		spt[0] = src[0] + hcell[0] * srcpts[i];
		for (j = 0; j < numsrcpts; ++j) {
			spt[1] = src[1] + hcell[1] * srcpts[j];
			for (l = 0; l < numsrcpts; ++l) {
				spt[2] = src[2] + hcell[2] * srcpts[l];
				val = grf (k, spt, obs);
				ans += srcwts[i] * srcwts[j] * srcwts[l] * val;
			}
		}
	}

	ans *= hcell[0] * hcell[1] * hcell[2];
	return ans;
}

/* N-point (per dimension) Duffy integration of the self term. The location src
 * is ignored because it is assumed to coincide with the origin. The location
 * obs specifies an observation (singularity) relative to the position of src
 * and is assumed to lie within the cell. */
cplx duffyint (real k, real *src, real *obs, real *cell, ifunc grf) {
	cplx ans = 0;
	int i;
	real ctr[3], dc[3], tmp;

	for (i = 0; i < 3; ++i) {
		/* Integrate over the pyramid along the +x axis. */
		duffcell (ctr, dc, obs, cell, 1);
		/* The observation is NULL because it is folded into the cell. */
		ans += srcint (k, ctr, NULL, dc, grf);

		/* Integrate over the pyramid along the -x axis. */
		duffcell (ctr, dc, obs, cell, -1);
		ans += srcint (k, ctr, NULL, dc, grf);

		/* Rotate the next axis into the x position. */
		rotate (obs);
		rotate (cell);
	}

	/* Completing the loop should restore the original order of obs and cell. */

	return ans;
}

void bldintrules (int nspts, int nrpts) {
	if (nspts > 0) {
		/* Allocate a source integration rule. */
		srcpts = malloc (2 * nspts * sizeof(real));
		srcwts = srcpts + nspts;
		gaussleg (srcpts, srcwts, nspts);
		numsrcpts = nspts;
	}
	if (nrpts > 0) {
		/* Allocate a separate receive integration rule. */
		rcvpts = malloc (2 * nrpts * sizeof(real));
		rcvwts = rcvpts + nrpts;
		gaussleg (rcvpts, rcvwts, nrpts);
		numrcvpts = nrpts;
	} else {
		/* If no receive rule was specified, reuse the source rule. */
		numrcvpts = numsrcpts;
		rcvpts = srcpts;
		rcvwts = srcwts;
	}
}

void delintrules () {
	numrcvpts = numsrcpts = 0;
	if (rcvpts && rcvpts != srcpts) free (rcvpts);
	if (srcpts) free (srcpts);

	srcpts = srcwts = rcvpts = rcvwts = NULL;
}

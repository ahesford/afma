#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "fsgreen.h"
#include "integrate.h"

void gaqd_ (int *, double *, double *, double *, double *, int *);

complex float radint (float k, float *cen, float *s, float *r,
		float *dc, int npts, double *pts, double *wts) {
	complex float ans = 0, val, csf;
	int i, j, l;
	float rv[3], sc, sr;

	sc = s[0] * cen[0] + s[1] * cen[1] + s[2] * cen[2];
	csf = cexp (I * k * sc);

	for (i = 0; i < npts; ++i) {
		rv[0] = r[0] + 0.5 * dc[0] * pts[i];
		for (j = 0; j < npts; ++j) {
			rv[1] = r[1] + 0.5 * dc[1] * pts[j];
			for (l = 0; l < npts; ++l) {
				rv[2] = r[2] + 0.5 * dc[2] * pts[l];
				sr = s[0] * rv[0] + s[1] * rv[1] + s[2] * rv[2];
				val = csf * cexp (-I * k * sr);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc[0] * dc[1] * dc[2] / 8;
	return ans;
}

complex float oneptint (float k, float *src, float *obs, float *dc) {
	complex float ans;
	float vol;

	vol = dc[0] * dc[1] * dc[2];

	ans = fsgreen (k, src, obs);
	ans *= vol;

	return vol;
}

/* N-point (per dimension) integration of both source and receiver. */
complex float srcint (float k, float *src, float *obs, float *dc,
		int nq, double *pts, double *wts) {
	complex float ans = 0, val;
	int i, j, l;
	float srcpt[3];

	for (i = 0; i < nq; ++i) {
		srcpt[0] = src[0] + 0.5 * dc[0] * pts[i];
		for (j = 0; j < nq; ++j) {
			srcpt[1] = src[1] + 0.5 * dc[1] * pts[j];
			for (l = 0; l < nq; ++l) {
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

/* High-accuracy evaluation of the self-interaction term. */
complex float selfinthigh (float k, float *dc) {
	complex float ans;
	float obspt[3] = { 0, 0, 0 };
	int nq = 16, ierr;
	double pts[16], wts[16];

	/* Get the quadrature points for the self-integration term. */
	gaqd_ (&nq, pts, wts, NULL, NULL, &ierr);
	for (ierr = 0; ierr < nq; ++ierr)
		pts[ierr] = cos (pts[ierr]);

	ans = srcint (k, obspt, obspt, dc, nq, pts, wts);
	ans *= k * k;

	return ans;
}

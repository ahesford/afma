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

/* N-point (per dimension) integration of the receiver. */
complex float rcvint (igrandf green, float k, float *src, float *cen, float *dc,
		int npts, double *pts, double *wts) {
	complex float ans = 0, val;
	int i, j, l;
	float obs[3];

	for (i = 0; i < npts; ++i) {
		obs[0] = cen[0] + 0.5 * dc[0] * pts[i];
		for (j = 0; j < npts; ++j) {
			obs[1] = cen[1] + 0.5 * dc[1] * pts[j];
			for (l = 0; l < npts; ++l) {
				obs[2] = cen[2] + 0.5 * dc[2] * pts[l];
				val = green (k, obs, src);
				ans += wts[i] * wts[j] * wts[l] * val;
			}
		}
	}

	ans *= dc[0] * dc[1] * dc[2] / 8;
	return ans;
}

/* One-point integration for the source, and four-point (per dimension)
 * integration of the receiver. Seems to be accurate enough... */
complex float fastint (float k, float *src, float *obs, float *dc,
		int npts, double *pts, double *wts) {
	complex float ans;

	/* Use four-point integration in the receiver box. */
	ans = rcvint (fsgreen, k, src, obs, dc, npts, pts, wts);
	ans *= dc[0] * dc[1] * dc[2];

	return ans;
}

complex float oneptint (float k, float *src, float *obs, float *dc) {
	complex float ans;
	float vol;

	vol = dc[0] * dc[1] * dc[2];

	ans = fsgreen (k, src, obs);
	ans *= vol * vol;

	return vol;
}

/* N-point (per dimension) integration of both source and receiver. */
complex float srcint (float k, float *src, float *obs, float *dc, int nsrc,
		double *spts, double *swts, int nrcv, double *rpts, double *rwts) {
	complex float ans = 0, val;
	int i, j, l;
	float srcpt[3];

	for (i = 0; i < nsrc; ++i) {
		srcpt[0] = src[0] + 0.5 * dc[0] * spts[i];
		for (j = 0; j < nsrc; ++j) {
			srcpt[1] = src[1] + 0.5 * dc[1] * spts[j];
			for (l = 0; l < nsrc; ++l) {
				srcpt[2] = src[2] + 0.5 * dc[2] * spts[l];
				val = rcvint (fsgreen, k, srcpt, obs, dc, nrcv, rpts, rwts);
				ans += swts[i] * swts[j] * swts[l] * val;
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
	float r, vol;

	vol = dc[0] * dc[1] * dc[2];
	r = cbrt (3 * vol / (4 * M_PI));
	ikr = I * k * r;

	ans = (1.0 - ikr) * cexp (ikr) - 1.0;
	ans *= vol;

	return ans;
}

/* High-accuracy evaluation of the self-interaction term. */
complex float selfinthigh (float k, float *dc) {
	complex float ans;
	float obspt[3] = { 0, 0, 0 };
	int nsrc = 16, nrcv = 14, ierr;
	double srcpts[16], srcwts[16], rcvpts[14], rcvwts[14];

	/* Get the quadrature points for the self-integration term. */
	gaqd_ (&nsrc, srcpts, srcwts, NULL, NULL, &ierr);
	for (ierr = 0; ierr < nsrc; ++ierr)
		srcpts[ierr] = cos (srcpts[ierr]);
	gaqd_ (&nrcv, rcvpts, rcvwts, NULL, NULL, &ierr);
	for (ierr = 0; ierr < nrcv; ++ierr)
		rcvpts[ierr] = cos (rcvpts[ierr]);

	ans = srcint (k, obspt, obspt, dc, nsrc, srcpts, srcwts, nrcv, rcvpts, rcvwts);
	ans *= k * k;

	return ans;
}

#include <complex.h>
#include <math.h>
#include <float.h>

#include <mpi.h>

#include "util.h"

/* The RMS error between a test vector and a reference. */
float mse (complex float *test, complex float *ref, long n, int nrm) {
	long i;
	float err[2] = {0.0, 0.0}, tmp;

	/* Calculate local contributions to the mean-squared error. */
	for (i = 0; i < n; ++i) {
		tmp = cabs(test[i] - ref[i]);
		err[0] += tmp * tmp;
		tmp = cabs(ref[i]);
		err[1] += tmp * tmp;
	}

	/* Sum the local contributions. */
	MPI_Allreduce (MPI_IN_PLACE, err, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	/* Return the normalized MSE if desired, otherwise just the difference. */
	if (nrm) return sqrt(err[0] / err[1]);
	return sqrt(err[0]);
}

/* Compute the sinc of the argument. */
float sinc (float x) {
	return ((fabs(x) < FLT_EPSILON) ? 1.0 : (sin (M_PI * x) / (M_PI * x)));
}

/* Compute the coordinates of the indexed far-field sample. */
int sampcoords (float *s, int i, float *th, int nt, int np) {
	int pi, ti;
	float st, theta, phi;

	/* The first sample is the south pole.
	 * No further computation is required. */
	if (i == 0) {
		s[0] = s[1] = 0.0;
		s[2] = -1.0;
		return 0;
	}

	/* Compute the sample indices in phi and theta. */
	pi = (i - 1) % np;
	ti = (i - 1) / np + 1;

	/* The sample is the north pole at the final theta index. */
	if (ti == nt - 1) {
		s[0] = s[1] = 0.0;
		s[2] = 1.0;
		return 0;
	}

	/* Compute the angular position from the sample indices. */
	phi = 2 * M_PI * pi / (float)np;
	theta = th[ti];

	/* Compute the cartesian position. */
	st = sin(theta);
	s[0] = cos(phi) * st;
	s[1] = sin(phi) * st;
	s[2] = cos(theta);

	return 0;
}

/* Compute the relative coordinates of a cell within a group of cells. */
int cellcoords (float *r, int l, int bpd, float dx) {
	int i, j, k;

	k = l % bpd;
	j = (l / bpd) % bpd;
	i = l / (bpd * bpd);

	r[0] = 0.5 * dx * (2.0 * i + 1.0 - (float)bpd);
	r[1] = 0.5 * dx * (2.0 * j + 1.0 - (float)bpd);
	r[2] = 0.5 * dx * (2.0 * k + 1.0 - (float)bpd);

	return 0;
}

/* Determine if the elemnt i is in the set s of length l. */
int inset (int i, int *s, int l) {
	int j;

	for (j = 0; j < l; ++j) 
		if (i == s[j]) return 1;

	return 0;
}

/* Find the index of the maximum value in set of length n, ignoring indices in
 * the set excl of length nex. */
int maxind (complex float *set, int n, int *excl, int nex) {
	float mv, cv;
	int mi, i;

	for (i = 0, mv = -1.0, mi = -1; i < n; ++i) {
		/* Skip the element if it is in the exclusion list. */
		if (inset(i, excl, nex)) continue;

		/* Determine if this value is the maximum. */
		cv = cabs(set[i]);
		if (cv > mv) {
			mi = i;
			mv = cv;
		}
	}

	return mi;
}
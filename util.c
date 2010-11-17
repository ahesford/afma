#include <complex.h>
#include <math.h>
#include <float.h>

#include <mpi.h>

#include "util.h"

/* Use the modified Gram-Schmidt process to compute (in place) the portion of
 * the n-dimensional vector v orthogonal to each of the nv vectors s. The
 * projection * of the vector onto each of the basis vectors is stored in the
 * length-nv array c. */
int cmgs (complex float *v, complex float *c, complex float *s,
		long n, int nv, int imgsit, float imgstol) {
	long i, j, k;
	complex float *sv, cv;

	for (i = 0, sv = s; i < nv; ++i, sv += n) {
		c[i] = 0;
		k = 0;

		do {
			cv = pardot (sv, v, n);
			c[i] += cv;
#pragma omp parallel for default(shared) private(j)
			for (j = 0; j < n; ++j) v[j] -= cv * sv[j];
		} while (cabs(cv / c[i]) > imgstol && ++k < imgsit);
		
	}

	return n;
}

/* Compute the inner product of the distributed vectors x and y of dimension n. */
complex float pardot (complex float *x, complex float *y, long n) {
	complex double dp = 0.0;
	long i;

#pragma omp parallel for default(shared) private(i) reduction(+: dp)
	for (i = 0; i < n; ++i) dp += conj(x[i]) * y[i];

	/* Add in the contributions from other processors. */
	MPI_Allreduce (MPI_IN_PLACE, &dp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return (complex float)dp;
}

float parnorm (complex float *x, long n) {
	double nrm = 0.0, nr, ni;
	long i;

#pragma omp parallel for default(shared) private(nr,ni,i) reduction(+: nrm)
	for (i = 0; i < n; ++i) {
		nr = creal(x[i]);
		nr *= nr;
		ni = cimag(x[i]);
		ni *= ni;
		nrm += nr + ni;
	}

	/* Sum over processors and reduce. */
	MPI_Allreduce (MPI_IN_PLACE, &nrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return (float)sqrt(nrm);
}

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

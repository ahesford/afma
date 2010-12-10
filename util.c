#include <mpi.h>

#include "precision.h"

#include "util.h"

/* Use the modified Gram-Schmidt process to compute (in place) the portion of
 * the n-dimensional vector v orthogonal to each of the nv vectors s. The
 * projection * of the vector onto each of the basis vectors is stored in the
 * length-nv array c. This is actually a selective reorthogonalization scheme
 * to attempt to correct rounding errors. */
int cmgs (cplx *v, cplx *c, cplx *s, long n, int nv) {
	long i, j;
	cplx *sv, cv;
	real vnrm, lcrit;

	/* Perform the first modified Gram Schmidt orthogonalization. */
	for (i = 0, sv = s, lcrit = 0.; i < nv; ++i, sv += n) {
		/* The projection of the vector onto the current basis. */
		c[i] = pardot (sv, v, n);

		/* Track the 1-norm of the projection column. */
		lcrit += cabs(c[i]);

#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < n; ++j) v[j] -= c[i] * sv[j];
	}

	/* Compute the norm of the vector. */
	c[nv] = parnorm (v, n);
	vnrm = creal(c[nv]);

	/* Reorthogonalize if necessary. */
	if (lcrit / vnrm > IMGS_L) {
		for (i = 0, sv = s; i < nv; ++i, sv += n)  {
			/* Re-project the vector onto the current basis. */
			cv = pardot (sv, v, n);

			/* Update the projection. */
			c[i] += cv;

			/* Remove the remaining parallel component. */
#pragma omp parallel for default(shared) private(j)
			for (j = 0; j < n; ++j) v[j] -= cv * sv[j];
		}

		/* Update the norm of the orthogonal vector. */
		c[nv] = parnorm (v, n);
		vnrm = creal(c[nv]);
	}

	/* Don't normalize if the norm is vanishing. */
	if (vnrm < REAL_EPSILON) return n;

	/* Finally, normalize the newly-created vector. */
#pragma omp parallel for default(shared) private(j)
	for (j = 0; j < n; ++j) v[j] /= vnrm;

	return n;
}

/* Compute the inner product of the distributed vectors x and y of dimension n. */
cplx pardot (cplx *x, cplx *y, long n) {
	/* Always compute the product in double precision to avoid rounding errors. */
	complex double dp = 0.0;
	long i;

#pragma omp parallel for default(shared) private(i) reduction(+: dp)
	for (i = 0; i < n; ++i) dp += conj(x[i]) * y[i];

	/* Add in the contributions from other processors. */
	MPI_Allreduce (MPI_IN_PLACE, &dp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return (cplx)dp;
}

real parnorm (cplx *x, long n) {
	/* Always compute the norm in double precision to avoid rounding errors. */
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

	return (real)sqrt(nrm);
}

/* The RMS error between a test vector and a reference. */
real mse (cplx *test, cplx *ref, long n, int nrm) {
	long i;
	real err[2] = {0.0, 0.0}, tmp;

	/* Calculate local contributions to the mean-squared error. */
	for (i = 0; i < n; ++i) {
		tmp = cabs(test[i] - ref[i]);
		err[0] += tmp * tmp;
		tmp = cabs(ref[i]);
		err[1] += tmp * tmp;
	}

	/* Sum the local contributions. */
	MPI_Allreduce (MPI_IN_PLACE, err, 2, MPIREAL, MPI_SUM, MPI_COMM_WORLD);

	/* Return the normalized MSE if desired, otherwise just the difference. */
	if (nrm) return sqrt(err[0] / err[1]);
	return sqrt(err[0]);
}

/* Compute the sinc of the argument. */
real sinc (real x) {
	return ((fabs(x) < REAL_EPSILON) ? 1.0 : (sin (M_PI * x) / (M_PI * x)));
}

/* Compute the coordinates of the indexed far-field sample. */
int sampcoords (real *s, int i, real *th, int nt, int np) {
	int pi, ti;
	real st, theta, phi;

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
	phi = 2 * M_PI * pi / (real)np;
	theta = th[ti];

	/* Compute the cartesian position. */
	st = sin(theta);
	s[0] = cos(phi) * st;
	s[1] = sin(phi) * st;
	s[2] = cos(theta);

	return 0;
}

/* Compute the relative coordinates of a cell within a group of cells. */
int cellcoords (real *r, int l, int bpd, real dx) {
	int idx[3];

	GRID(idx, l, bpd, bpd);

	r[0] = 0.5 * dx * (2.0 * idx[0] + 1.0 - (real)bpd);
	r[1] = 0.5 * dx * (2.0 * idx[1] + 1.0 - (real)bpd);
	r[2] = 0.5 * dx * (2.0 * idx[2] + 1.0 - (real)bpd);

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
int maxind (cplx *set, int n, int *excl, int nex) {
	real mv, cv;
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

/* Calculate the Legendre polynomial of order m and its derivative at a point t. */
static int legendre (real *p, real *dp, real t, int m) {
	real p0 = 1.0, p1 = t;
	int k;

	/* Handle low orders explicitly. */
	if (m < 1) {
		*p = 1.0;
		*dp = 0.0;
		return m;
	} else if (m < 2) {
		*p = t;
		*dp = 1.0;
		return m;
	}

	/* Otherwise, calculate the function values recursively. */
	for (k = 1; k < m; ++k) {
		*p = ((2.0 * k + 1.0) * t * p1 - k * p0) / (1.0 + k);
		p0 = p1; p1 = *p;
	}

	*dp = m * (p0 - t * p1) / (1.0 - t * t);

	return m;
}

/* Compute Gaussian quadrature nodes and weights. */
int gaussleg (real *nodes, real *weights, int m) {
	int i, j, nroots = (m + 1) / 2;
	real t, p, dp, dt;
	const real tol = REAL_EPSILON;
	const int maxit = 100;

	for (i = 0; i < nroots; ++i) {
		/* Use the Chebyshev roots as an initial guess. */
		t = cos (M_PI * (i + 0.75) / (m + 0.5));
		for (j = 0; j < maxit; ++j) {
			/* Compute the value of the Legendre polynomial. */
			legendre (&p, &dp, t, m);
			/* Perform a Newton-Raphson update. */
			dt = -p / dp; t += dt;

			/* Break if convergence detected. */
			if (fabs(dt) < tol) break;
		}

		/* Update the nodes and weights. */
		nodes[i] = t;
		nodes[m - i - 1] = -t;
		weights[i] = 2.0 / (1.0 - t * t) / (dp * dp);
		weights[m - i - 1] = weights[i];
	}

	return 0;
}

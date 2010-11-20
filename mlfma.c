#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include <mpi.h>

/* Pull in the CBLAS header. */
#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#else
#ifdef _FREEBSD
#include <cblas.h>
#else
#include <gsl_cblas.h>
#endif /* _FREEBSD */
#endif /* _MACOSX */

#include "ScaleME.h"

#include "io.h"
#include "mlfma.h"
#include "direct.h"
#include "util.h"
#include "fsgreen.h"
#include "integrate.h"

fmadesc fmaconf;

static int farmatrow (complex float *, float, float *, float, int);
static int farmatcol (complex float *, float, float *, float *, int, int);
static int acabuild (complex float **, float, float, float *, int, int, float, int);
static int fullbuild (complex float **, float, float *, int, int, float, int);
static int recaca (complex float *, complex float *, int, int, int, float);

/* Computes the far-field pattern for the specified group with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. The list of "basis
 * functions" and center are ignored, since the calling program ensures each
 * FMM basis function is actually a single group of gridded elements with the
 * same affine grid. The center of the finest-level FMM group coincides with
 * the single FMM basis function contained therein. */
void farpattern (int nbs, int *bsl, void *vcrt, void *vpat, float *cen, int sgn) {
	complex float fact, beta = 1.0, *crt = (complex float *)vcrt,
		*pat = *((complex float **)vpat);

	if (sgn >= 0) {
		/* Scalar factors for the matrix multiplication. */
		fact = fmaconf.k0 * fmaconf.cellvol;

		/* Don't add in the existing output vector since it should
		 * have been zeroed anyway. */
		beta = 0.0;

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasNoTrans, fmaconf.nsamp,
				fmaconf.bspboxvol, &fact, fmaconf.radpats,
				fmaconf.nsamp, crt, 1, &beta, pat, 1);
	} else {
		/* Distribute the far-field patterns to the basis functions. */
		/* Scalar factors for the matrix multiplication. */
		fact = I * fmaconf.k0 * fmaconf.k0 / (4 * M_PI);

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasConjTrans, fmaconf.nsamp,
				fmaconf.bspboxvol, &fact, fmaconf.radpats,
				fmaconf.nsamp, pat, 1, &beta, crt, 1);
	}
}

/* Computes the far-field pattern for the specified group with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. The list of "basis
 * functions" and center are ignored, since the calling program ensures each
 * FMM basis function is actually a single group of gridded elements with the
 * same affine grid. The center of the finest-level FMM group coincides with
 * the single FMM basis function contained therein. ACA is used to approximate
 * the matrix for efficient computations. */
void acafarpattern (int nbs, int *bsl, void *vcrt, void *vpat, float *cen, int sgn) {
	complex float fact, beta = 1.0, *crt = (complex float *)vcrt,
		*pat = *((complex float **)vpat), *work, *u, *v;

	u = fmaconf.radpats;
	v = fmaconf.radpats + fmaconf.acarank * fmaconf.nsamp;

	work = malloc (fmaconf.acarank * sizeof(complex float));

	if (sgn >= 0) {
		beta = 0.0;
		fact = 1.0;

		cblas_cgemv (CblasColMajor, CblasConjTrans, fmaconf.bspboxvol,
				fmaconf.acarank, &fact, v, fmaconf.bspboxvol,
				crt, 1, &beta, work, 1);

		/* Scalar factors for the matrix multiplication. */
		fact = fmaconf.k0 * fmaconf.cellvol;

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasNoTrans, fmaconf.nsamp,
				fmaconf.acarank, &fact, u, fmaconf.nsamp,
				work, 1, &beta, pat, 1);
	} else {
		beta = 0.0;
		fact = 1.0;

		cblas_cgemv (CblasColMajor, CblasConjTrans, fmaconf.nsamp,
				fmaconf.acarank, &fact, u, fmaconf.nsamp,
				pat, 1, &beta, work, 1);

		/* Distribute the far-field patterns to the basis functions. */
		/* Scalar factors for the matrix multiplication. */
		fact = I * fmaconf.k0 * fmaconf.k0 / (4 * M_PI);
		beta = 1.0;

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasNoTrans, fmaconf.bspboxvol,
				fmaconf.acarank, &fact, v, fmaconf.bspboxvol,
				work, 1, &beta, crt, 1);
	}

	free (work);
}

/* Build a column of the far-field signature for a point a distance rmc
 * from the center of the parent box. */
static int farmatcol (complex float *col, float k, float *rmc,
		float *thetas, int ntheta, int nphi) {
	int nsamp = nphi * (ntheta - 2) + 2;

#pragma omp parallel default(shared)
{
	float s[3];
	int i;
#pragma omp for
	for (i = 0; i < nsamp; ++i) {
		/* Compute the Cartesian coordinates of the far-field sample. */
		sampcoords (s, i, thetas, ntheta, nphi);
		/* Compute the far-field sample. */
		col[i] = fsplane (k, rmc, s);
	}
}

	return 0;
}

/* Build a row of the far-field signature for a fixed field sample s. */
static int farmatrow (complex float *row, float k, float *s, float dx, int bpd) {
	int bptot = bpd * bpd * bpd;

#pragma omp parallel default(shared)
{
	int l;
	float dist[3];

#pragma omp for
	for (l = 0; l < bptot; ++l) {
		/* The distance from the basis to the box center. */
		cellcoords (dist, l, bpd, dx);

		/* The value of the far-field sample of the cell. */
		row[l] = fsplane (k, dist, s);
	}
}
	return 0;
}

/* Construct an ACA approximation to the far-field signature matrix. */
static int acabuild (complex float **mats, float k0, float tol, float *thetas,
		int ntheta, int nphi, float dx, int bpd) {
	int maxrank, *irow, *icol, i, j, k,
	    nsamp = (ntheta - 2) * nphi + 2, nelt = bpd * bpd * bpd;
	complex float *u, *v, *uptr, *vptr, *row, *col, dpr, dpc;
	float dist[3], s[3], err = 0, tolsq;

	/* Compute the square of the tolerance and decrease by two orders
	 * of magnitude for recompression. */
	tolsq = tol * tol * 1e-2;

	maxrank = MIN(nelt, nsamp);

	/* Allocate the row and column index arrays. The extra element should
	 * only be stored (but never recalled) in the limiting case when the
	 * rank cannot be reduced. */
	irow = malloc((2 * maxrank + 1) * sizeof(int));
	icol = irow + maxrank + 1;

	/* Allocate the workspace for the row and column matrices. */
	u = malloc(maxrank * (nelt + nsamp) * sizeof(complex float));
	v = u + nsamp * maxrank;

	/* Start with the first row of the matrix. */
	irow[0] = 0;

	for (i = 0, row = v, col = u; i < maxrank; ++i, row += nelt, col += nsamp) {
		/* Find the coordinate of the observer element. */
		sampcoords (s, irow[i], thetas, ntheta, nphi);

		/* Fill the row for the selected observer element. */
		farmatrow (row, k0, s, dx, bpd);

		/* Subtract existing contributions from earlier ranks. */
#pragma omp parallel for default(shared) private(vptr, uptr, j, k)
		for (j = 0; j < nelt; ++j) {
			for (k = 0; k < i; ++k, uptr += nsamp, vptr += nelt) {
				uptr = u + k * nsamp;
				vptr = v + k * nelt;
				row[j] -= vptr[j] * uptr[irow[i]];
			}
		}

		/* Find the index and value of the maximum value. */
		icol[i] = maxind (row, nelt, icol, i);
		dpr = row[icol[i]];

		/* Scale the row. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) row[j] /= dpr;

		/* Find the coordinates of the source cell. */
		cellcoords (dist, icol[i], bpd, dx);
	
		/* Construct the radiation pattern of the source cell. */
		farmatcol (col, k0, dist, thetas, ntheta, nphi);

		/* Subtract existing contributions from earlier ranks. */
#pragma omp parallel for default(shared) private(vptr, uptr, j, k)
		for (j = 0; j < nsamp; ++j) {
			for (k = 0; k < i; ++k) {
				vptr = v + k * nelt;
				uptr = u + k * nsamp;
				col[j] -= vptr[icol[i]] * uptr[j];
			}
		}

		/* Update the error approximation. */
		for (k = 0, uptr = u, vptr = v; k < i; ++k) {
			cblas_cdotu_sub (nelt, vptr, 1, row, 1, &dpr);
			cblas_cdotu_sub (nsamp, uptr, 1, col, 1, &dpc);
			uptr += nsamp;
			vptr += nelt;

			err += 2.0 * cabs(dpr) * cabs(dpc);
		}

		cblas_cdotc_sub (nelt, row, 1, row, 1, &dpr);
		cblas_cdotc_sub (nsamp, col, 1, col, 1, &dpc);

		err += creal(dpr) * creal(dpc);

		if (creal(dpr) * creal(dpc) <= tolsq * err) break;

		/* Find the next row index. */
		irow[i + 1] = maxind (col, nsamp, irow, i + 1);
	}

	maxrank = i;

	/* Conjugate the matrix v. */
	for (i = 0, k = maxrank * nelt; i < k; ++i) v[i] = conj(v[i]);

	/* Recompress (in place) the matrices with the actual tolerance. */
	maxrank = recaca (u, v, nsamp, nelt, maxrank, tol);

	/* Allocate the final matrix storage. */
	*mats = malloc (maxrank * (nelt + nsamp) * sizeof(complex float));
	/* Copy the colum matrix in first, then the row matrix. */
	memcpy (*mats, u, maxrank * nsamp * sizeof(complex float));
	memcpy (*mats + maxrank * nsamp, v, maxrank * nelt * sizeof(complex float));

	/* Free the work arrays. */
	free (irow);
	free (u);

	return maxrank;
}

/* Recompress an ACA approximation using truncated singular values. */
static int recaca (complex float *u, complex float *v, int m, int n, int k, float tol) {
	int rank, lwork, info, i, j, l;
	float *rwork, *ss;
	complex float *qu, *qv, *rp, *work, *tu, *tv, *us, *vs, alpha = 1.0, beta = 0.0;

	/* Allocate space for the orthogonal matrices and the outer product. */
	qu = malloc (k * (m + n + 3 * k) * sizeof(complex float));
	qv = qu + k * m;
	rp = qv + k * n;
	us = rp + k * k;
	vs = us + k * k;

	/* Allocate space for the reflector coefficients. */
	tu = malloc (2 * k * sizeof(complex float));
	tv = tu + k;

	/* Allocate space for the singular values and real work. */
	ss = malloc (6 * k * sizeof(float));
	rwork = ss + k;

	/* Figure out the maximum size of the work array. */
	lwork = -1;
	cgeqrf_ (&m, &k, qu, &m, tu, qu, &lwork, &info);
	cgeqrf_ (&n, &k, qv, &n, tv, qu + 1, &lwork, &info);
	cungqr_ (&m, &k, &k, qu, &m, tu, qu + 2, &lwork, &info);
	cungqr_ (&n, &k, &k, qv, &n, tv, qu + 3, &lwork, &info);
	cgesvd_ ("S", "S", &k, &k, rp, &k, ss, us, &k, vs, &k, qu + 4, &lwork, rwork, &info);

	/* Find the maximum workspace size. */
	for (i = 0; i < 5; ++i) lwork = MAX(lwork, (int)creal(qu[i]));

	work = malloc(lwork * sizeof(complex float));

	/* Copy the input matrices into the workspaces. */
	memcpy (qu, u, (m * k) * sizeof(complex float));
	memcpy (qv, v, (n * k) * sizeof(complex float));

	/* Compute the QR factorizations of the row and column matrices. */
	cgeqrf_ (&m, &k, qu, &m, tu, work, &lwork, &info);
	cgeqrf_ (&n, &k, qv, &n, tv, work, &lwork, &info);

	/* Compute the outer product of the triangular matrices. */
	for (i = 0; i < k; ++i) {
		for (j = 0; j < k; ++j) {
			rp[i + j * k] = 0.0;
			for (l = MAX(i,j); l < k; ++l)
				rp[i + j * k] += qu[i + l * m] * conj(qv[j + l * n]);
		}
	}

	/* Compute the SVD of the triangular matrix product. */	
	cgesvd_ ("S", "S", &k, &k, rp, &k, ss, us, &k, vs, &k, work, &lwork, rwork, &info);

	/* Find the truncated rank. */
	for (rank = 0; rank < k; ++rank)
		if (fabs(ss[rank] / ss[0]) < tol) break;

	/* Multiply the singular values by the column matrix. */
	for (i = 0; i < rank; ++i)
		for (j = 0; j < k; ++j)
			vs[i + k * j] *= ss[i];

	/* Expand the orthogonal matrices. */
	cungqr_ (&m, &k, &k, qu, &m, tu, work, &lwork, &info);
	cungqr_ (&n, &k, &k, qv, &n, tv, work, &lwork, &info);

	/* Compute the new row and column matrices. */
	cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, m, rank, k,
			&alpha, qu, m, us, k, &beta, u, m);
	cblas_cgemm (CblasColMajor, CblasNoTrans, CblasConjTrans, n, rank, k,
			&alpha, qv, n, vs, k, &beta, v, n);

	free (qu);
	free (tu);
	free (ss);
	free (work);

	return rank;
}

static int svdbuild (complex float **mats, float k0, float tol, float *thetas,
		int ntheta, int nphi, float dx, int bpd) {
	int nsamp = (ntheta - 2) * nphi + 2, nelt = bpd * bpd * bpd,
	    lwork, info, mindim, rank, l, i;
	complex float *col, *a, *u, *vt, *work, *ptr;
	float dist[3], *s, *rwork;

	/* The minimum dimension. */
	mindim = MIN(nsamp, nelt);

	/* Allocate the full matrix. */
	a = malloc(nsamp * nelt * sizeof(complex float));

	/* Allocate the singular vector matrices. */
	u = malloc(mindim * (nsamp + nelt) * sizeof(complex float));
	vt = u + mindim * nsamp;

	/* Allocate the real workspace and the singular value array. */
	s = malloc(6 * mindim * sizeof(float));
	rwork = s + mindim;

	/* Perform a workspace query. */
	lwork = -1;
	cgesvd_ ("S", "S", &nsamp, &nelt, a, &nsamp, s, u, &nsamp, vt, &mindim,
			a, &lwork, rwork, &info);

	/* Allocate the workspace. */
	lwork = creal(a[0]);
	work = malloc(lwork * sizeof(complex float));

	/* Loop through all columns (source grid elements) and build the matrix. */
	for (l = 0, col = a; l < nelt; ++l, col += nsamp) {
		/* The relative position of the source grid element. */
		cellcoords (dist, l, bpd, dx);

		/* Build the corresponding matrix column. */
		farmatcol (col, k0, dist, thetas, ntheta, nphi);
	}

	/* Perform an SVD on the matrix. */
	cgesvd_ ("S", "S", &nsamp, &nelt, a, &nsamp, s, u, &nsamp, vt, &mindim,
			work, &lwork, rwork, &info);

	/* Find the first rank below the desired tolerance. */
	for (rank = 0; rank < mindim; ++rank)
		if (fabs(s[rank] / s[0]) < tol) break;

	*mats = malloc (rank * (nelt + nsamp) * sizeof(complex float));

	/* Multiply the left singular vector matrix by the singular values. */
	for (l = 0, col = u, ptr = *mats; l < rank; ++l)
		for (i = 0; i < nsamp; ++i, ++col, ++ptr)
			*ptr = s[l] * (*col);

	/* Tranpose and conjugate the right singular vector matrix. */
	for (l = 0, ptr = *mats + rank * nsamp; l < rank; ++l)
		for (i = 0; i < nelt; ++i)
			ptr[i + l * nelt] = conj(vt[l + i * mindim]);

	free (a);
	free (u);
	free (s);
	free (work);

	return rank;
}

static int fullbuild (complex float **mats, float k0, float *thetas,
		int ntheta, int nphi, float dx, int bpd) {
	int nsamp = (ntheta - 2) * nphi + 2, nelt = bpd * bpd * bpd, l;
	complex float *col;
	float dist[3];

	/* Allocate the full far-field matrix. */
	*mats = malloc(nsamp * nelt * sizeof(complex float));

	/* Loop through all columns (source grid elements). */
	for (l = 0, col = *mats; l < nelt; ++l, col += nsamp) {
		/* The relative position of the source grid element. */
		cellcoords (dist, l, bpd, dx);

		/* Build the corresponding matrix column. */
		farmatcol (col, k0, dist, thetas, ntheta, nphi);
	}

	return nelt * nsamp;
}

/* Precomputes the near interactions for redundant calculations and sets up
 * the wave vector directions to be used for fast calculation of far-field patterns. */
int fmmprecalc (float acatol, int useaca) {
	float *thetas;
	int ntheta, nphi, rank, i;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Get the finest level parameters. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, NULL, NULL);
	/* Allocate the theta array. */
	thetas = malloc (ntheta * sizeof(float));
	/* Populate the theta array. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, thetas, NULL);
	
	/* Convert the theta samples from the cosines of the angles to the angles. */
	for (i = 0; i < ntheta; ++i) thetas[i] = acos(thetas[i]);


	if (!useaca) {
		/* Build the direct far-field matrices. */
		fmaconf.acarank = 0;
		i = fullbuild (&(fmaconf.radpats), fmaconf.k0, thetas,
				ntheta, nphi, fmaconf.cell, fmaconf.bspbox);
		fprintf (stderr, "Rank %d: Far-field matrix element count: %d\n", rank, i);
	} else {
		/* Use ACA if the tolerance is positive, otherwise use SVD. */
		if (acatol > 0)
			fmaconf.acarank = acabuild (&(fmaconf.radpats),
					fmaconf.k0, acatol, thetas, ntheta,
					nphi, fmaconf.cell, fmaconf.bspbox);
		else
			fmaconf.acarank = svdbuild (&(fmaconf.radpats),
					fmaconf.k0, -acatol, thetas, ntheta,
					nphi, fmaconf.cell, fmaconf.bspbox);
		
		fprintf (stderr, "Rank %d: Far-field matrix rank: %d\n", rank, fmaconf.acarank);
	}

	free (thetas);

	return fmaconf.acarank;
}

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (int useaca) {
	int error;
	float len, cen[3];
	
	/* The problem and tree are both three-dimensional. */
	ScaleME_setDimen (3);
	ScaleME_setTreeType (3);

	/* The fields are scalar-valued. */
	ScaleME_setFields (1);

	/* The wave number is real-valued. */
	ScaleME_setWaveNumber (fmaconf.k0);

	/* Set some MLFMA parameters. */
	ScaleME_setNumBasis (fmaconf.gnumbases);
	ScaleME_setMaxLevel (fmaconf.maxlev);
	ScaleME_setPrecision (fmaconf.precision);
	ScaleME_setMAC (fmaconf.numbuffer);
	ScaleME_setInterpOrder (fmaconf.interpord);

	ScaleME_setTopComputeLevel (fmaconf.toplev);

	ScaleME_setRHSDataWidth (fmaconf.bspboxvol);

	if (fmaconf.fo2iterm > 0) {
		ScaleME_useFastO2IGhosts (1);
		ScaleME_selectFastO2I (fmaconf.fo2iterm, fmaconf.fo2iord, fmaconf.fo2iosr);
	}

	/* Set the root box length. */
	len = (1 <<  fmaconf.maxlev) * fmaconf.bspbox * fmaconf.cell;
	/* Position the root box properly. */
	cen[0] = fmaconf.min[0] + 0.5 * len;
	cen[1] = fmaconf.min[1] + 0.5 * len;
	cen[2] = fmaconf.min[2] + 0.5 * len;
	ScaleME_setRootBox (len, cen);
	
	/*  let all processes start the initialisation together */
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* open the std files */
	MPFMA_stdout = stdout;
	MPFMA_stderr = stderr;

	/* Use the external near-field interactions. */
	ScaleME_setBlockDirInterFunc (blockinteract);
	if (useaca) ScaleME_useExternFarField (acafarpattern);
	else ScaleME_useExternFarField (farpattern);

	/* Finish the setup with the external interactions. */
	error = ScaleME_initSetUp (MPI_COMM_WORLD, NULL, NULL, NULL, bscenter);

	if (error) {
		fprintf(stdout, "ERROR: ScaleME pre-init failed.\n");
		ScaleME_finalizeParHostFMA();
		return -1;
	}

	return 0;
}

int ScaleME_postconf (void) {
	if (ScaleME_completeSetUp()) {
		fprintf(stdout, "ERROR: ScaleME setup routine failed.\n");
		goto error_handle;
	}

	if (ScaleME_initParHostDataStructs()) {
		fprintf(stdout, "ERROR: ScaleME parallel host init failed.\n");
		goto error_handle;
	}

	return 0;

error_handle:
	ScaleME_finalizeParHostFMA();
	return -1;
}

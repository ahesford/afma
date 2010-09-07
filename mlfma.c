#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <mpi.h>

#include <omp.h>

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

#include "mlfma.h"
#include "direct.h"
#include "util.h"
#include "fsgreen.h"
#include "integrate.h"

fmadesc fmaconf;

static int farmatrow (complex float *, float, float *, float, int);
static int farmatcol (complex float *, float, float *, float *, int, int);

/* Computes the far-field pattern for the specified group with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. The list of "basis
 * functions" and center are ignored, since the calling program ensures each
 * FMM basis function is actually a single group of gridded elements with the
 * same affine grid. The center of the finest-level FMM group coincides with
 * the single FMM basis function contained therein. */
void farpattern (int nbs, int *bsl, void *vcrt, void *vpat, float *cen, int sgn) {
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
int acabuild (complex float **mats, float k0, float tol, float *thetas,
		int ntheta, int nphi, float dx, int bpd) {
	int maxrank, *irow, *icol, i, j, k,
	    nsamp = (ntheta - 2) * nphi + 2, nelt = bpd * bpd * bpd;
	complex float *u, *v, *uptr, *vptr, *row, *col, dpr, dpc;
	float dist[3], s[3], err = 0;

	tol *= tol;

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

		if (creal(dpr) * creal(dpc) <= tol * err) break;

		/* Find the next row index. */
		irow[i + 1] = maxind (col, nsamp, irow, i + 1);
	}

	maxrank = i;

	/* Allocate the final matrix storage. */
	*mats = malloc (maxrank * (nelt * nsamp) * sizeof(complex float));
	/* Copy the colum matrix in first, then the row matrix. */
	memcpy (*mats, u, maxrank * nsamp * sizeof(complex float));
	/* The row matrix should be conjugated for ease of application. */
	for (i = 0, k = maxrank * nelt, j = maxrank * nsamp; i < k; ++i)
		(*mats)[j + i] = conj(v[i]);

	/* Free the work arrays. */
	free (irow);
	free (u);

	return maxrank;
}

/* Precomputes the near interactions for redundant calculations and sets up
 * the wave vector directions to be used for fast calculation of far-field patterns. */
int fmmprecalc () {
	float *thetas;
	int ntheta, nphi, rank, i;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Get the finest level parameters. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, NULL, NULL);
	/* Allocate the theta array. */
	thetas = malloc (ntheta * sizeof(float));
	/* Populate the theta array. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, thetas, NULL);

	/* Convert the theta samples to angular values rather than cosine angles. */
	for (i = 0; i < ntheta; ++i) thetas[i] = acos(thetas[i]);

	/* Build the ACA approximation to far-field matrices. */
	fmaconf.acarank = acabuild (&(fmaconf.radpats), fmaconf.k0, 1e-6, thetas,
			ntheta, nphi, fmaconf.cell, fmaconf.bspbox);

	free (thetas);

	fprintf (stderr, "Rank %d: Radiation pattern matrix rank: %d\n", rank, fmaconf.acarank);

	return fmaconf.acarank;
}

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (void) {
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
	ScaleME_useExternFarField (farpattern);

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

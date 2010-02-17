#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <mpi.h>

#include <omp.h>

#ifdef _MACOSX
#include <Accelerate/Accelerate.h>
#endif /* _MACOSX */

#ifdef _FREEBSD
#include <cblas.h>
#endif /* _FREEBSD */

#include "ScaleME.h"

#include "mlfma.h"
#include "direct.h"
#include "util.h"
#include "fsgreen.h"
#include "integrate.h"

fmadesc fmaconf;

/* Computes the far-field pattern for the specified basis with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. */
void farpattern (int nbs, int *bsl, void *vcrt, void *vpat, float *cen, int sgn) {
	complex float fact, beta, *buf, *crt = (complex float *)vcrt,
		*pat = *((complex float **)vpat);
	int l, idx[3];

	/* The buffer that will be used by this thread for this routine. */
	buf = calloc (fmaconf.bspboxvol, sizeof(complex float));

	if (sgn >= 0) {
		/* Scalar factors for the matrix multiplication. */
		fact = fmaconf.k0;
		beta = 1.0;

		/* Compute the far-field pattern for the basis functions. */
		for (l = 0; l < nbs; ++l) {
			/* Find the basis grid position. */
			bsindex (bsl[l], idx);

			/* Shift to a local grid position. */
			idx[0] %= fmaconf.bspbox;
			idx[1] %= fmaconf.bspbox;
			idx[2] %= fmaconf.bspbox;
			
			/* Augment the current vector. */
			buf[SQIDX(fmaconf.bspbox,idx[0],idx[1],idx[2])] = crt[l];
		}

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasNoTrans, fmaconf.nsamp,
				fmaconf.bspboxvol, &fact, fmaconf.radpats, 
				fmaconf.nsamp, buf, 1, &beta, pat, 1);
	} else {
		/* Distribute the far-field patterns to the basis functions. */
		/* Scalar factors for the matrix multiplication. */
		fact = I * fmaconf.k0 * fmaconf.k0 / (4 * M_PI);
		beta = 0.0;

		/* Perform the matrix-vector product. */
		cblas_cgemv (CblasColMajor, CblasConjTrans, fmaconf.nsamp,
				fmaconf.bspboxvol, &fact, fmaconf.radpats,
				fmaconf.nsamp, pat, 1, &beta, buf, 1);

		for (l = 0; l < nbs; ++l) {
			/* Find the basis grid position. */
			bsindex (bsl[l], idx);
			
			/* Shift to a local grid position. */
			idx[0] %= fmaconf.bspbox;
			idx[1] %= fmaconf.bspbox;
			idx[2] %= fmaconf.bspbox;
			
			/* Augment the current vector. */
			crt[l] += buf[SQIDX(fmaconf.bspbox,idx[0],idx[1],idx[2])];
		}
	}

	free (buf);
}

/* Precompute the exponential radiation pattern for a point a distance rmc 
 * from the center of the parent box. */
int buildradpat (complex float *pat, float k, float *rmc,
		float *thetas, int ntheta, int nphi) {
	int i, j, nthsc = ntheta - 1;
	float s[3], dphi = 2 * M_PI / nphi, sn, phi;

	/* South pole first. */
	s[0] = s[1] = 0.0;
	s[2] = -1.0;
	*(pat++) = srcint (k, rmc, s, fmaconf.cell, fsplane);

	for (i = 1; i < nthsc; ++i) {
		s[2] = thetas[i];
		sn = sin (acos (thetas[i]));

		for (j = 0, phi = 0; j < nphi; ++j, phi += dphi) {
			s[0] = sn * cos (phi);
			s[1] = sn * sin (phi);
			*(pat++) = srcint (k, rmc, s, fmaconf.cell, fsplane);
		}
	}

	/* North pole last. */
	s[0] = s[1] = 0.0;
	s[2] = 1.0;
	*pat = srcint (k, rmc, s, fmaconf.cell, fsplane);

	return ntheta;
}

/* Precomputes the near interactions for redundant calculations and sets up
 * the wave vector directions to be used for fast calculation of far-field patterns. */
int fmmprecalc () {
	float *thetas, clen;
	int ntheta, nphi, rank;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Get the finest level parameters. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, NULL, NULL);
	/* Allocate the theta array. */
	thetas = malloc (ntheta * sizeof(float));
	/* Populate the theta array. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, thetas, NULL);

	/* Allocate storage for the radiation patterns. */
	fmaconf.radpats = malloc (fmaconf.bspboxvol * fmaconf.nsamp * sizeof(complex float));

	fprintf (stderr, "Rank %d: Radiation pattern buffer size: %ld bytes\n",
			rank, fmaconf.bspboxvol * fmaconf.nsamp * sizeof(complex float));

	/* Calculate the box center. */
	clen = 0.5 * (float)fmaconf.bspbox;

#pragma omp parallel default(shared)
{
	int i, j, k, l;
	complex float *pptr;
	float dist[3];

#pragma omp for
	for (l = 0; l < fmaconf.bspboxvol; ++l) {
		/* The pointer to the relevant pattern. */
		pptr = fmaconf.radpats + l * fmaconf.nsamp;

		/* The basis index with respect to the parent box. */
		k = l % fmaconf.bspbox;
		j = (l / fmaconf.bspbox) % fmaconf.bspbox;
		i = l / (fmaconf.bspbox * fmaconf.bspbox);

		/* The distance from the basis to the box center. */
		dist[0] = ((float)i + 0.5 - clen) * fmaconf.cell;
		dist[1] = ((float)j + 0.5 - clen) * fmaconf.cell;
		dist[2] = ((float)k + 0.5 - clen) * fmaconf.cell;

		/* Construct the radiation pattern. */
		buildradpat (pptr, fmaconf.k0, dist, thetas, ntheta, nphi);
	}
}

	free (thetas);

	return fmaconf.nsamp;
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

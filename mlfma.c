#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include <mpi.h>

#include "ScaleME.h"

#include "fsgreen.h"
#include "integrate.h"
#include "itsolver.h"
#include "mlfma.h"

#define UNSGN(x) ((x) > 0 ? (x) : -(x))

fmadesc fmaconf;

/* Computes the far-field pattern for the specified basis with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. */
void farpattern (void *vcrt, void *vpat, int gi, int sgn, float *cen) {
	float *s, sdcr, rv[3], fact;
	complex float *crt = (complex float *)vcrt, *pat = *((complex float **)vpat);
	int i;

	/* Make sure the sign has unity magnitude. */
	sgn = ((sgn >= 0) ? 1 : -1);

	/* Find the basis center. */
	bscenter (gi, rv);

	/* Compute the far-field pattern for the basis function. */
	if (sgn > 0) {
		fact = fmaconf.k0 * fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];
		for (i = 0, s = fmaconf.kvecs; i < fmaconf.nsamp; ++i, s+= 3) {
			sdcr = s[0] * (cen[0] - rv[0]) + s[1] * (cen[1] - rv[1])
				+ s[2] * (cen[2] - rv[2]);
			
			pat[i] += *crt * fact * cexp (I * fmaconf.k0 * sdcr);
		}
	} else {
		fact = fmaconf.k0 * fmaconf.k0 / (4 * M_PI);
		for (i = 0, s = fmaconf.kvecs; i < fmaconf.nsamp; ++i, s += 3) {
			sdcr = s[0] * (rv[0] - cen[0]) + s[1] * (rv[1] - cen[1])
				+ s[2] * (rv[2] - cen[2]);

			*crt += I * pat[i] * fact * cexp (I * fmaconf.k0 * sdcr);
		}
	}
}

/* Precomputes the near interactions for redundant calculations and sets up
 * the wave vector directions to be used for fast calculation of far-field patterns. */
int fmmprecalc () {
	float lbox[3], dist[3], zero[3] = { 0, 0, 0 }, *kptr, *thetas, phi, dphi, sn;
	int i, j, k, idx, totel, nb1, ntheta, nphi;

	/* Get the finest level parameters. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, NULL, NULL);

	/* Allocate the theta array. */
	thetas = malloc (ntheta * sizeof(float));
	fmaconf.kvecs = malloc (3 * fmaconf.nsamp * sizeof(float));

	/* Populate the theta array. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, thetas, NULL);

	/* The step in phi. */
	dphi = 2 * M_PI / nphi;

	/* The south pole comes first. */
	kptr = fmaconf.kvecs;
	kptr[0] = kptr[1] = 0.0;
	kptr[2] = -1.0;
	kptr += 3;

	/* The intermediate values are next. */
	for (i = 1; i < ntheta - 1; ++i) {
		sn = sin (acos (thetas[i]));
		for (j = 0, phi = 0; j < nphi; ++j, phi += dphi, kptr += 3) {
			kptr[0] = sn * cos (phi);
			kptr[1] = sn * sin (phi);
			kptr[2] = thetas[i];
		}
	}

	/* Now the north pole. */
	kptr[0] = kptr[1] = 0.0;
	kptr[2] = 1.0;

	/* Find the smallest box size. */
	ScaleME_getSmallestBoxSize (lbox);

	nb1 = fmaconf.numbuffer + 1;

	/* Find the maximum number of near neighbors in each dimension. */
	fmaconf.nbors[0] = (int) ((nb1 * lbox[0]) / fmaconf.cell[0]) + 1;
	fmaconf.nbors[1] = (int) ((nb1 * lbox[1]) / fmaconf.cell[1]) + 1;
	fmaconf.nbors[2] = (int) ((nb1 * lbox[2]) / fmaconf.cell[2]) + 1;

	totel = fmaconf.nbors[0] * fmaconf.nbors[1] * fmaconf.nbors[2];

	fmaconf.gridints = malloc (totel * sizeof (complex float));

	/* Compute the interactions. */
	for (i = 0, idx = 0; i < fmaconf.nbors[0]; ++i) {
		dist[0] = i * fmaconf.cell[0];
		for (j = 0; j < fmaconf.nbors[1]; ++j) {
			dist[1] = j * fmaconf.cell[1];
			for (k = 0; k < fmaconf.nbors[2]; ++k, ++idx) {
				dist[2] = k * fmaconf.cell[2];
				
				fmaconf.gridints[idx] = srcint (fmaconf.k0, zero, dist, fmaconf.cell);
				fmaconf.gridints[idx] *= fmaconf.k0 * fmaconf.k0;
			}
		}
	}

	/* The self term uses the analytic approximation. */
	fmaconf.gridints[0] = selfint (fmaconf.k0, fmaconf.cell);

	free (thetas);

	return totel;
}

/* Evaluate at a group of observers the fields due to a group of sources. */
void blockinteract (int nsrc, int nobs, int *srclist,
		int *obslist, void *vsrc, void *vobs) {
	int i, j, *src, *srcptr, obs[3], dist[3], idx;
	complex float *srcfld, *csrc, *cobs;

	/* Convert the void types to complex floats. */
	csrc = (complex float *)vsrc;
	cobs = (complex float *)vobs;

	src = malloc (3 * nsrc * sizeof(int));

	/* Precompute the source indices. */
	for (i = 0; i < nsrc; ++i) bsindex (srclist[i], src + 3 * i);

	for (i = 0; i < nobs; ++i, ++cobs) {
		/* Find the observer index. */
		bsindex (obslist[i], obs);

		/* Point to the starting source. */
		srcptr = src;
		srcfld = csrc;

		for (j = 0; j < nsrc; ++j, srcptr += 3) {
			/* Compute the grid distance between source and observer. */
			dist[0] = srcptr[0] - obs[0];
			dist[1] = srcptr[1] - obs[1];
			dist[2] = srcptr[2] - obs[2];
			
			/* Find the linear index for the interaction. */
			idx = UNSGN(dist[2]) + UNSGN(dist[1]) * fmaconf.nbors[2]
				+ UNSGN(dist[0]) * fmaconf.nbors[2] * fmaconf.nbors[1];

			/* Augment the observer field with this contribution.
			 * Also shift the pointers for the next iteration. */
			*cobs += *(srcfld++) * fmaconf.gridints[idx];
		}
	}

	free (src);

	return;
}

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (void) {
	int error;
	
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

	if (fmaconf.fo2iterm > 0)
		ScaleME_selectFastO2I (fmaconf.fo2iterm, fmaconf.fo2iord, fmaconf.fo2iosr);

	if (fmaconf.smallbox > 0)
		ScaleME_setSmallestBoxSize(fmaconf.smallbox);
	
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

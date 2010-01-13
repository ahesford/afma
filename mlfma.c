#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <mpi.h>

#include "ScaleME.h"

#include "fsgreen.h"
#include "integrate.h"
#include "itsolver.h"
#include "mlfma.h"

#define UNSGN(x) ((x) > 0 ? (x) : -(x))
#define IDX(nx,i,j,k) ((k) + (nx) * ((j) + (nx) * (i)))

fmadesc fmaconf;

/* Computes the far-field pattern for the specified basis with the specified
 * center, and stores the output in a provided vector. sgn is positive for
 * radiation pattern and negative for receiving pattern. */
void farpattern (void *vcrt, void *vpat, int gi, int sgn, float *cen) {
	float rv[3], fact, shift;
	complex float *crt = (complex float *)vcrt,
		*pat = *((complex float **)vpat), *ppat;
	int i, j, k;

	shift = 0.5 * ((float)fmaconf.bspbox - 1.0);

	/* Find the basis center. */
	bscenter (gi, rv);

	/* Find the position relative to the parent. */
	i = round((rv[0] - cen[0]) / fmaconf.cell + shift);
	j = round((rv[1] - cen[1]) / fmaconf.cell + shift);
	k = round((rv[2] - cen[2]) / fmaconf.cell + shift);

	/* Pull out the relevant precomputed pattern. */
	k += j * fmaconf.bspbox + i * fmaconf.bspbox * fmaconf.bspbox;
	ppat = fmaconf.radpats[k];

	/* Compute the far-field pattern for the basis function. */
	if (sgn >= 0) {
		fact = fmaconf.k0 * fmaconf.cellvol;
		for (i = 0; i < fmaconf.nsamp; ++i)
			pat[i] += *crt * fact * conj(ppat[i]);
	} else {
		fact = fmaconf.k0 * fmaconf.k0 / (4 * M_PI);
		for (i = 0; i < fmaconf.nsamp; ++i)
			*crt += I * pat[i] * fact * ppat[i];
	}
}

/* Precompute the exponential radiation pattern for a point a distance rmc 
 * from the center of the parent box. */
int buildradpat (complex float *pat, float k, float *rmc,
		float *thetas, int ntheta, int nphi) {
	int i, j, nthsc = ntheta - 1;
	float s[3], sdr, dphi = 2 * M_PI / nphi, sn, phi;

	/* South pole first. */
	sdr = -rmc[2];
	*(pat++) = cexp (I * k * sdr);

	for (i = 1; i < nthsc; ++i) {
		s[2] = thetas[i];
		sn = sin (acos (thetas[i]));

		for (j = 0, phi = 0; j < nphi; ++j, phi += dphi) {
			s[0] = sn * cos (phi);
			s[1] = sn * sin (phi);
			sdr = s[0] * rmc[0] + s[1] * rmc[1] + s[2] * rmc[2];
			*(pat++) = cexp (I * k * sdr);
		}
	}

	/* North pole last. */
	sdr = rmc[2];
	*pat = cexp (I * k * sdr);

	return ntheta;
}

/* Precomputes the near interactions for redundant calculations and sets up
 * the wave vector directions to be used for fast calculation of far-field patterns. */
int fmmprecalc () {
	float dist[3], zero[3] = { 0, 0, 0 }, *thetas, clen;
	int i, j, k, l, ip, jp, kp, totel, ntheta, nphi, nchild;
	complex float **pptr;

	/* The number of children in a finest-level box. */
	nchild = fmaconf.bspbox * fmaconf.bspbox * fmaconf.bspbox;

	/* Get the finest level parameters. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, NULL, NULL);
	/* Allocate the theta array. */
	thetas = malloc (ntheta * sizeof(float));
	/* Populate the theta array. */
	ScaleME_getFinestLevelParams (&(fmaconf.nsamp), &ntheta, &nphi, thetas, NULL);

	/* Allocate storage for the radiation patterns. */
	fmaconf.radpats = malloc (nchild * sizeof(complex float));
	fmaconf.radpats[0] = malloc (nchild * fmaconf.nsamp * sizeof(complex float));
	for (i = 1; i < nchild; ++i)
		fmaconf.radpats[i] = fmaconf.radpats[i - 1] + fmaconf.nsamp;

	/* Calculate the box center. */
	clen = 0.5 * (float)fmaconf.bspbox;

	/* Compute radiation patterns for children bases. */
	for (pptr = fmaconf.radpats, i = 0; i < fmaconf.bspbox; ++i) {
		dist[0] = ((float)i + 0.5 - clen) * fmaconf.cell;
		for (j = 0; j < fmaconf.bspbox; ++j) {
			dist[1] = ((float)j + 0.5 - clen) * fmaconf.cell;
			for (k = 0; k < fmaconf.bspbox; ++k, ++pptr) {
				dist[2] = ((float)k + 0.5 - clen) * fmaconf.cell;
				buildradpat (*pptr, fmaconf.k0, dist,
						thetas, ntheta, nphi);
			}
		}
	}

	/* Find the maximum number of near neighbors in each dimension. */
	fmaconf.nbors = (2 * fmaconf.numbuffer + 1) * fmaconf.bspbox;
	fmaconf.nborsex = 2 * fmaconf.nbors;

	/* Build the expanded grid. */
	totel = fmaconf.nborsex * fmaconf.nborsex * fmaconf.nborsex;
	fmaconf.gridints = fftw_malloc (totel * sizeof (complex double));

	fmaconf.fplan = fftw_plan_dft_3d (fmaconf.nborsex, fmaconf.nborsex,
			fmaconf.nborsex, fmaconf.gridints, fmaconf.gridints,
			FFTW_FORWARD, FFTW_MEASURE);
	fmaconf.bplan = fftw_plan_dft_3d (fmaconf.nborsex, fmaconf.nborsex,
			fmaconf.nborsex, fmaconf.gridints, fmaconf.gridints,
			FFTW_BACKWARD, FFTW_MEASURE);

	/* Compute the interactions. */
	for (l = 0; l < totel; ++l) {
		k = l % fmaconf.nborsex;
		j = (l / fmaconf.nborsex) % fmaconf.nborsex;
		i = l / (fmaconf.nborsex * fmaconf.nborsex);

		ip = (i < fmaconf.nbors) ? i : (fmaconf.nborsex - i);
		jp = (j < fmaconf.nbors) ? j : (fmaconf.nborsex - j);
		kp = (k < fmaconf.nbors) ? k : (fmaconf.nborsex - k);

		dist[0] = (float)ip * fmaconf.cell;
		dist[1] = (float)jp * fmaconf.cell;
		dist[2] = (float)kp * fmaconf.cell;
		
		fmaconf.gridints[l] = srcint (fmaconf.k0, zero, dist, fmaconf.cell);
		fmaconf.gridints[l] *= fmaconf.k0 * fmaconf.k0 / totel;
	}

	/* The self term uses the analytic approximation. */
	fmaconf.gridints[0] = selfint (fmaconf.k0, fmaconf.cell) / totel;

	/* Perform the Fourier transform of the Green's function. */
	fftw_execute (fmaconf.fplan);

	free (thetas);

	return fmaconf.nbors;
}

/* Evaluate at a group of observers the fields due to a group of sources. */
void blockinteract (int nsrc, int nobs, int *srclist,
		int *obslist, void *vsrc, void *vobs, float *bc) {
	int i, j, k, l, idx, totel;
	complex float *csrc = (complex float *)vsrc, *cobs = (complex float *)vobs;
	complex double *buf;
	float cc[3], shift;

	/* Allocate and clear the buffer array. */
	totel = fmaconf.nborsex * fmaconf.nborsex * fmaconf.nborsex;
	buf = fftw_malloc (totel * sizeof(complex double));
	memset (buf, 0, totel * sizeof(complex double));

	/* An offset to compute the grid position. */
	shift = 0.5 * (float)(fmaconf.nbors - 1);

	/* Populate the local grid. */
	for (l = 0; l < nsrc; ++l) {
		bscenter (srclist[l], cc);
		i = round((cc[0] - bc[0]) / fmaconf.cell + shift);
		j = round((cc[1] - bc[1]) / fmaconf.cell + shift);
		k = round((cc[2] - bc[2]) / fmaconf.cell + shift);
		idx = IDX(fmaconf.nborsex, i, j, k);

		buf[idx] = csrc[l];
	}

	/* Transform the local grid in place. */
	fftw_execute_dft (fmaconf.fplan, buf, buf);

	/* The convolution is now a multiplication. */
	for (l = 0; l < totel; ++l) buf[l] *= fmaconf.gridints[l];

	/* Inverse transform the grid in place. */
	fftw_execute_dft (fmaconf.bplan, buf, buf);

	/* Augment with output with the local convolution. */
	for (l = 0; l < nobs; ++l) {
		bscenter (obslist[l], cc);
		i = round((cc[0] - bc[0]) / fmaconf.cell + shift);
		j = round((cc[1] - bc[1]) / fmaconf.cell + shift);
		k = round((cc[2] - bc[2]) / fmaconf.cell + shift);
		idx = IDX(fmaconf.nborsex, i, j, k);

		cobs[l] += buf[idx];
	}

	/* Free the buffer array. */
	fftw_free (buf);

	return;
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

	if (fmaconf.fo2iterm > 0)
		ScaleME_selectFastO2I (fmaconf.fo2iterm, fmaconf.fo2iord, fmaconf.fo2iosr);

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

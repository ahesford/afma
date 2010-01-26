#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>

#include <omp.h>

#include "ScaleME.h"

#include "direct.h"
#include "util.h"
#include "fsgreen.h"
#include "mlfma.h"
#include "integrate.h"

dirdesc dircache;

complex float *dirbuf;
int totbpnbr;

/* Initialize the direct-interaction cache structure. */
int mkdircache () {
	int nbs, *bslist, l, idx[3];

	/* Get the complete list of locally required basis functions. */
	ScaleME_getLocallyReqBasis (&nbs, &bslist);

	/* Get the index of the first basis function. */
	bsindex (bslist[0], idx);
	/* Get the box index for this basis function. */
	idx[0] /= fmaconf.bspbox;
	idx[1] /= fmaconf.bspbox;
	idx[2] /= fmaconf.bspbox;

	/* Set the box bounds for the local RHS cache. */
	dircache.boxmin[0] = dircache.boxmax[0] = idx[0];
	dircache.boxmin[1] = dircache.boxmax[1] = idx[1];
	dircache.boxmin[2] = dircache.boxmax[2] = idx[2];

	/* Now loop through and establish the maxmimum and minimum box bounds. */
	for (l = 1; l < nbs; ++l) {
		/* Get the index of the first basis function. */
		bsindex (bslist[l], idx);
		
		/* Get the box index for this basis function. */
		idx[0] /= fmaconf.bspbox;
		idx[1] /= fmaconf.bspbox;
		idx[2] /= fmaconf.bspbox;

		/* Update the maximum and minimum box bounds. */
		dircache.boxmin[0] = MIN(dircache.boxmin[0], idx[0]);
		dircache.boxmin[1] = MIN(dircache.boxmin[1], idx[1]);
		dircache.boxmin[2] = MIN(dircache.boxmin[2], idx[2]);
		dircache.boxmax[0] = MAX(dircache.boxmax[0], idx[0]);
		dircache.boxmax[1] = MAX(dircache.boxmax[1], idx[1]);
		dircache.boxmax[2] = MAX(dircache.boxmax[2], idx[2]);
	}

	/* Count the boxes in the box cache. */
	dircache.nbox[0] = dircache.boxmax[0] - dircache.boxmin[0];
	dircache.nbox[1] = dircache.boxmax[1] - dircache.boxmin[1];
	dircache.nbox[2] = dircache.boxmax[2] - dircache.boxmin[2];
	dircache.nboxprod = dircache.nbox[0] * dircache.nbox[1] * dircache.nbox[2];

	/* The total number of direct interaction elements in the cache. */
	dircache.totelts = 8 * dircache.nboxprod * fmaconf.bspboxvol;

	/* Allocate the cache arrays. */
	dircache.boxfill = malloc (dircache.nboxprod * sizeof(int));
	dircache.boxrhs = malloc (dircache.totelts * sizeof(complex float));

	/* The locally-required basis list is no longer necessary. */
	free (bslist);

	return dircache.nboxprod;
}

/* Clear the cache arrays for a new iteration. */
void clrdircache () {
	memset (dircache.boxfill, 0, dircache.nboxprod * sizeof(int));
	memset (dircache.boxrhs, 0,  dircache.totelts * sizeof(complex float));
}

/* Free the allocated memory in the direct-interaction cache structure. */
void freedircache () {
	free (dircache.boxfill);
	free (dircache.boxrhs);

	dircache.nbox[2] = dircache.nbox[1] = dircache.nbox[0] =
		dircache.nboxprod = dircache.totelts = 0;
}

/* Precompute some values for the direct interactions. */
int dirprecalc () {
	/* The number of neighbor boxes per dimension. */
	dircache.nbors = dircache.nborsvol = 2 * fmaconf.numbuffer + 1;
	/* The number of near-neighbor boxes total. */
	dircache.nborsvol *= dircache.nborsvol * dircache.nborsvol;

	/* The FFT size. */
	dircache.nfft[0] = dircache.nfft[1] = dircache.nfft[2] = 2 * fmaconf.bspbox;
	dircache.nfftprod = dircache.nfft[0] * dircache.nfft[1] * dircache.nfft[2];

	/* Build the expanded grid. */
	totbpnbr = dircache.nfftprod * dircache.nborsvol;
	dircache.gridints = fftwf_malloc (totbpnbr * (1 + omp_get_max_threads()) * sizeof(complex float));
	dirbuf = dircache.gridints + totbpnbr;

	/* The forward FFT plan transforms all boxes in one pass. */
	dircache.fplan = fftwf_plan_many_dft (3, dircache.nfft, dircache.nborsvol,
			dircache.gridints, NULL, 1, dircache.nfftprod, dircache.gridints,
			NULL, 1, dircache.nfftprod, FFTW_FORWARD, FFTW_MEASURE);
	/* The inverse FFT plan only transforms a single box. */
	dircache.bplan = fftwf_plan_dft_3d (dircache.nfft[0], dircache.nfft[1], dircache.nfft[2],
			dircache.gridints, dircache.gridints, FFTW_BACKWARD, FFTW_MEASURE);

#pragma omp parallel default(shared)
{
	int off[3], l, i, j, k;
	complex float *grf;

#pragma omp for
	for (l = 0; l < dircache.nborsvol; ++l) {
		k = l % dircache.nbors;
		j = (l / dircache.nbors) % dircache.nbors;
		i = l / (dircache.nbors * dircache.nbors);

		grf = dircache.gridints + l * dircache.nfftprod;

		off[0] = (i - fmaconf.numbuffer) * fmaconf.bspbox;
		off[1] = (j - fmaconf.numbuffer) * fmaconf.bspbox;
		off[2] = (k - fmaconf.numbuffer) * fmaconf.bspbox;

		/* Build the Green's function grid for this local box. */
		greengrid (grf, fmaconf.bspbox, dircache.nfft[0], fmaconf.k0, fmaconf.cell, off);
	}
}

	/* Perform the Fourier transform of the Green's function. */
	fftwf_execute (dircache.fplan);

	return dircache.nfftprod;
}

/* Evaluate at a group of observers the fields due to a group of sources. */
void blockinteract (int nsrc, int nobs, int *srclist,
		int *obslist, void *vsrc, void *vobs, float *bc) {
	int l, bsoff[3], idx[3], i, j, k;
	complex float *csrc = (complex float *)vsrc, *cobs = (complex float *)vobs;
	complex float *buf, *bptr, *cbox;
	float fbox[3];

	/* Allocate and clear the buffer array. */
	buf = dirbuf + omp_get_thread_num() * totbpnbr;
	memset (buf, 0, totbpnbr * sizeof(complex float));

	/* Calculate the box position. */
	fbox[0] = (bc[0] - fmaconf.min[0]) / (fmaconf.cell * fmaconf.bspbox);
	fbox[1] = (bc[1] - fmaconf.min[1]) / (fmaconf.cell * fmaconf.bspbox);
	fbox[2] = (bc[2] - fmaconf.min[2]) / (fmaconf.cell * fmaconf.bspbox);

	/* Now calculate an offset for the basis indices. */
	bsoff[0] = ((int)(fbox[0]) - fmaconf.numbuffer) * fmaconf.bspbox;
	bsoff[1] = ((int)(fbox[1]) - fmaconf.numbuffer) * fmaconf.bspbox;
	bsoff[2] = ((int)(fbox[2]) - fmaconf.numbuffer) * fmaconf.bspbox;

	/* Populate the local grid. */
	for (l = 0; l < nsrc; ++l) {
		bsindex (srclist[l], idx);

		/* Convert the global grid position to a local position. */
		idx[0] -= bsoff[0];
		idx[1] -= bsoff[1];
		idx[2] -= bsoff[2];

		/* Find the local box number. */
		i = idx[0] / fmaconf.bspbox;
		j = idx[1] / fmaconf.bspbox;
		k = idx[2] / fmaconf.bspbox;

		/* Point to the local buffer for this box. */
		bptr = buf + dircache.nfftprod * SQIDX(dircache.nbors,i,j,k);

		/* Find the position in the local box. */
		i = idx[0] % fmaconf.bspbox;
		j = idx[1] % fmaconf.bspbox;
		k = idx[2] % fmaconf.bspbox;

		bptr[SQIDX(dircache.nfft[0],i,j,k)] = csrc[l];
	}

	/* Transform the local grids in place. */
	fftwf_execute_dft (dircache.fplan, buf, buf);

	/* The convolutions are now multiplications. */
	for (l = 0; l < totbpnbr; ++l) buf[l] *= dircache.gridints[l];

	/* Point to the center box. */
	cbox = buf + dircache.nfftprod * fmaconf.numbuffer * 
		(1 + dircache.nbors * (1 + dircache.nbors));
	for (l = 0, bptr = buf; l < dircache.nborsvol; ++l, bptr += dircache.nfftprod) {
		if (bptr == cbox) continue;
		for (i = 0; i < dircache.nfftprod; ++i) cbox[i] += bptr[i];
	}

	/* Inverse transform the grid in place. */
	fftwf_execute_dft (dircache.bplan, cbox, cbox);

	/* Augment with output with the local convolution. */
	for (l = 0; l < nobs; ++l) {
		bsindex (obslist[l], idx);

		/* Convert the global grid position to a local position. */
		idx[0] %= fmaconf.bspbox;
		idx[1] %= fmaconf.bspbox;
		idx[2] %= fmaconf.bspbox;

		cobs[l] += cbox[SQIDX(dircache.nfft[0],idx[0],idx[1],idx[2])];
	}

	return;
}

/* Build the extended Green's function on an expanded cubic grid. */
int greengrid (complex float *grf, int m, int mex, float k0, float cell, int *off) {
	int i, j, k, ip, jp, kp;
	float dist[3], zero[3] = {0., 0., 0.}, scale;

	/* The scale of the integral equation solution. */
	scale = k0 * k0 / (float)(mex * mex * mex);

	/* Compute the interactions. */
	for (i = 0; i < mex; ++i) {
		ip = (i < m) ? i : (i - mex);
		dist[0] = (float)(ip - off[0]) * fmaconf.cell;
		for (j = 0; j < mex; ++j) {
			jp = (j < m) ? j : (j - mex);
			dist[1] = (float)(jp - off[1]) * fmaconf.cell;
			for (k = 0; k < mex; ++k) {
				kp = (k < m) ? k : (k - mex);
				dist[2] = (float)(kp - off[2]) * fmaconf.cell;

				/* Handle the self term specially. */
				if (kp == off[2] && jp == off[1] && ip == off[0])
					*(grf++) = selfint (k0, cell) / (mex * mex * mex);
				else *(grf++) = scale * srcint (k0, zero, dist, cell);
			}
		}
	}

	return mex;
}

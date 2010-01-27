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

/* Buffers for the RHS cache, the Green's functions, and a workspace. */
complex float *dirbuf, **boxrhs, *gridints;
int ncbox[3], cboxmin[3], *boxfill;
int nbors, nfftprod, nfft[3];
fftwf_plan fplan, bplan;
omp_lock_t **boxlocks;

/* Initialize the direct-interaction cache structure. */
int mkdircache () {
	int nbs, *bslist, l, i, idx[3], boxmax[3], nebox, ncboxprod;

	/* Get the complete list of locally required basis functions. */
	ScaleME_getLocallyReqBasis (&nbs, &bslist);

	/* The total number of nonempty boxes in the locally-required list. */
	nebox = GEDIV(nbs,fmaconf.bspboxvol);

	/* Get the index of the first basis function. */
	bsindex (bslist[0], idx);
	/* Get the box index for this basis function. */
	idx[0] /= fmaconf.bspbox;
	idx[1] /= fmaconf.bspbox;
	idx[2] /= fmaconf.bspbox;

	/* Set the box bounds for the local RHS cache. */
	cboxmin[0] = boxmax[0] = idx[0];
	cboxmin[1] = boxmax[1] = idx[1];
	cboxmin[2] = boxmax[2] = idx[2];

	/* Now loop through and establish the maxmimum and minimum box bounds. */
	for (l = 1; l < nbs; ++l) {
		/* Get the index of the first basis function. */
		bsindex (bslist[l], idx);
		
		/* Get the box index for this basis function. */
		idx[0] /= fmaconf.bspbox;
		idx[1] /= fmaconf.bspbox;
		idx[2] /= fmaconf.bspbox;

		/* Update the maximum and minimum box bounds. */
		cboxmin[0] = MIN(cboxmin[0], idx[0]);
		cboxmin[1] = MIN(cboxmin[1], idx[1]);
		cboxmin[2] = MIN(cboxmin[2], idx[2]);
		boxmax[0] = MAX(boxmax[0], idx[0]);
		boxmax[1] = MAX(boxmax[1], idx[1]);
		boxmax[2] = MAX(boxmax[2], idx[2]);
	}

	/* Count the boxes in the box cache. */
	ncbox[0] = boxmax[0] - cboxmin[0] + 1;
	ncbox[1] = boxmax[1] - cboxmin[1] + 1;
	ncbox[2] = boxmax[2] - cboxmin[2] + 1;
	ncboxprod = ncbox[0] * ncbox[1] * ncbox[2];

	/* Allocate the cache arrays. */
	boxfill = calloc (ncboxprod, sizeof(int));
	boxrhs = calloc (ncboxprod, sizeof(complex float *));
	boxrhs[0] = malloc (nebox * nfftprod * sizeof(complex float));

	/* Allocate box locks. */
	boxlocks = calloc (ncboxprod, sizeof(omp_lock_t *));
	boxlocks[0] = malloc (nebox * sizeof(omp_lock_t));

	/* Now loop through all bases and set up the box cache pointers. */
	for (l = 0, nebox = 0; l < nbs; ++l) {
		/* Get the index of the first basis function. */
		bsindex (bslist[l], idx);
		
		/* Get the box index for this basis function. */
		idx[0] = (idx[0] / fmaconf.bspbox) - cboxmin[0];
		idx[1] = (idx[1] / fmaconf.bspbox) - cboxmin[1];
		idx[2] = (idx[2] / fmaconf.bspbox) - cboxmin[2];

		/* Find the linear box index. */
		i = IDX(ncbox[1], ncbox[2], idx[0], idx[1], idx[2]);

		/* Skip this basis function if the box has already been counted. */
		if (boxfill[i]) continue;

		/* Mark this box as counted. */
		boxfill[i] = 1;

		/* Set the next box pointer and initialize the lock. */
		boxrhs[i] = boxrhs[0] + nebox * nfftprod;
		boxlocks[i] = boxlocks[0] + nebox;
		omp_init_lock (boxlocks[i]);
		++nebox;
	}

	/* The locally-required basis list is no longer necessary. */
	free (bslist);

	return ncboxprod;
}

/* Clear the cache arrays for a new iteration. */
void clrdircache () {
	int ncboxprod = ncbox[0] * ncbox[1] * ncbox[2];
	memset (boxfill, 0, ncboxprod * sizeof(int));
}

/* Free the allocated memory in the direct-interaction cache structure. */
void freedircache () {
	free (boxfill);
	free (boxrhs[0]);
	free (boxrhs);
	free (gridints);
	free (dirbuf);
	free (boxlocks[0]);
	free (boxlocks);
}

/* Check the cache for the given RHS. If it exists, return the pre-cached copy.
 * Otherwise, cache the provided copy and take the DFT before returning a pointer
 * to the cache bin. */
complex float *cacheboxrhs (complex float *rhs, int *bslist, int nbs) {
	int idx[3], l, i;
	complex float *bptr;

	/* Get the index for the first basis in the box. */
	bsindex (bslist[0], idx);

	/* The local box index. */
	idx[0] = (idx[0] / fmaconf.bspbox) - cboxmin[0];
	idx[1] = (idx[1] / fmaconf.bspbox) - cboxmin[1];
	idx[2] = (idx[2] / fmaconf.bspbox) - cboxmin[2];

	/* Get the index in the cache. */
	l = IDX(ncbox[1], ncbox[2], idx[0], idx[1], idx[2]);
	/* Point to the storage in the cache. */
	bptr = boxrhs[l];

	/* The cache check and fill operation must be thread safe. */
	omp_set_lock (boxlocks[l]);

	/* Cache miss. Fill the box. */
	if (!(boxfill[l])) {
		/* Clear the cache storage. */
		memset (bptr, 0, nfftprod * sizeof(complex float)); 
		
		/* Populate the local grid. */
		for (i = 0; i < nbs; ++i) {
			/* Find the basis index. */
			bsindex (bslist[i], idx);
			
			/* Find the position in the local box. */
			idx[0] %= fmaconf.bspbox;
			idx[1] %= fmaconf.bspbox;
			idx[2] %= fmaconf.bspbox;
			
			/* Fill the expanded FFT grid. */
			bptr[SQIDX(nfft[0],idx[0],idx[1],idx[2])] = rhs[i];
		}
		
		/* Transform the cached RHS. */
		fftwf_execute_dft (fplan, bptr, bptr);
		
		/* Mark the cache spot as full. */
		boxfill[l] = 1;
	}

	/* Free the lock to allow other threads to access the cache block. */
	omp_unset_lock (boxlocks[l]);

	/* Return the RHS. */
	return bptr;
}

/* Precompute some values for the direct interactions. */
int dirprecalc () {
	int totbpnbr, nborsvol;

	/* The number of neighbor boxes per dimension. */
	nbors = nborsvol = 2 * fmaconf.numbuffer + 1;
	/* The number of near-neighbor boxes total. */
	nborsvol *= nborsvol * nborsvol;

	/* The FFT size. */
	nfft[0] = nfft[1] = nfft[2] = 2 * fmaconf.bspbox;
	nfftprod = nfft[0] * nfft[1] * nfft[2];

	/* Build the expanded grid. */
	totbpnbr = nfftprod * nborsvol;
	gridints = fftwf_malloc (totbpnbr * sizeof(complex float));
	dirbuf = fftwf_malloc (nfftprod * omp_get_max_threads() * sizeof(complex float));

	/* The forward FFT plan transforms all boxes in one pass. */
	fplan = fftwf_plan_dft_3d (nfft[0], nfft[1], nfft[2],
			gridints, gridints, FFTW_FORWARD, FFTW_MEASURE);
	/* The inverse FFT plan only transforms a single box. */
	bplan = fftwf_plan_dft_3d (nfft[0], nfft[1], nfft[2],
			gridints, gridints, FFTW_BACKWARD, FFTW_MEASURE);

#pragma omp parallel default(shared)
{
	int off[3], l, i, j, k;
	complex float *grf;

#pragma omp for
	for (l = 0; l < nborsvol; ++l) {
		k = l % nbors;
		j = (l / nbors) % nbors;
		i = l / (nbors * nbors);

		grf = gridints + l * nfftprod;

		off[0] = (i - fmaconf.numbuffer) * fmaconf.bspbox;
		off[1] = (j - fmaconf.numbuffer) * fmaconf.bspbox;
		off[2] = (k - fmaconf.numbuffer) * fmaconf.bspbox;

		/* Build the Green's function grid for this local box. */
		greengrid (grf, fmaconf.bspbox, nfft[0], fmaconf.k0, fmaconf.cell, off);

		/* Fourier transform the Green's function. */
		fftwf_execute_dft (fplan, grf, grf);
	}
}

	/* Allocate the local cache structure. */
	mkdircache ();

	return nfftprod;
}

/* Evaluate at a group of observers the fields due to a group of sources. */
void blockinteract (int nsrc, int nobs, int *srclist,
		int *obslist, void *vsrc, void *vobs) {
	int i, l, boxoff[3], idx[3];
	complex float *csrc = (complex float *)vsrc, *cobs = (complex float *)vobs;
	complex float *buf, *gptr, *bptr;

	/* Clear the local output buffer. */
	buf = dirbuf + omp_get_thread_num() * nfftprod;
	memset (buf, 0, nfftprod * sizeof(complex float));

	/* Find the index for the first basis in the target box. */
	bsindex (obslist[0], boxoff);

	/* Find the minimum box index for near interactions. */
	boxoff[0] = (boxoff[0] / fmaconf.bspbox) - fmaconf.numbuffer;
	boxoff[1] = (boxoff[1] / fmaconf.bspbox) - fmaconf.numbuffer;
	boxoff[2] = (boxoff[2] / fmaconf.bspbox) - fmaconf.numbuffer;

	/* Populate the local grid. */
	for (l = 0; l < nsrc; l += fmaconf.bspboxvol) {
		/* Get the cached RHS for the source box in question. */
		bptr = cacheboxrhs (csrc + l, srclist + l, fmaconf.bspboxvol);

		/* Get the index fo the first basis in the source box. */
		bsindex (srclist[l], idx);

		/* Find the local box number. */
		idx[0] = (idx[0] / fmaconf.bspbox) - boxoff[0];
		idx[1] = (idx[1] / fmaconf.bspbox) - boxoff[1];
		idx[2] = (idx[2] / fmaconf.bspbox) - boxoff[2];

		/* Point to the Green's function for this box. */
		gptr = gridints + nfftprod * SQIDX(nbors,idx[0],idx[1],idx[2]);

		/* Convolve the source field with the Green's function and
		 * augment the field at the target. */
		for (i = 0; i < nfftprod; ++i) buf[i] += gptr[i] * bptr[i];
	}

	/* Inverse transform the grid in place. */
	fftwf_execute_dft (bplan, buf, buf);

	/* Augment with output with the local convolution. */
	for (l = 0; l < nobs; ++l) {
		bsindex (obslist[l], idx);

		/* Convert the global grid position to a local position. */
		idx[0] %= fmaconf.bspbox;
		idx[1] %= fmaconf.bspbox;
		idx[2] %= fmaconf.bspbox;

		/* Augment the RHS. */
		cobs[l] += buf[SQIDX(nfft[0],idx[0],idx[1],idx[2])];
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

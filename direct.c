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

typedef struct {
	int index[3], fill;
	omp_lock_t lock;
	complex float *rhs;
} boxdesc;

/* Buffers for the RHS cache, the Green's functions, and a workspace. */
static boxdesc *boxlist;
static complex float *gridints, *rhsbuf;
static int nbors, nfftprod, nfft[3], nebox;
static fftwf_plan fplan, bplan;

/* Compare two box indices for sorting and searching. */
int idxcomp (const void *vl, const void *vr) {
	int *bl = (int *)vl, *br = (int *)vr;

	if (bl[0] != br[0]) return bl[0] - br[0];
	if (bl[1] != br[1]) return bl[1] - br[1];
	return bl[2] - br[2];
}

/* Compare two box indices within a box descriptor structure. */
int boxcomp (const void *vl, const void *vr) {
	boxdesc *bl = (boxdesc *)vl, *br = (boxdesc *)vr;

	if (bl->index[0] != br->index[0]) return bl->index[0] - br->index[0];
	if (bl->index[1] != br->index[1]) return bl->index[1] - br->index[1];
	return bl->index[2] - br->index[2];
}

/* Initialize the direct-interaction cache structure. */
int mkdircache () {
	complex float *rhsptr;
	int nbs, *bslist, i, *idx, *boxidx, rank;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Get the complete list of locally required basis functions. */
	ScaleME_getLocallyReqBasis (&nbs, &bslist);

	/* Allocate the box index array. */
	boxidx = malloc (3 * nbs * sizeof(int));

	/* Build the array of box indices. */
	for (i = 0, idx = boxidx; i < nbs; ++i, idx += 3) {
		bsindex (bslist[i], idx);
		idx[0] /= fmaconf.bspbox;
		idx[1] /= fmaconf.bspbox;
		idx[2] /= fmaconf.bspbox;
	}

	/* The basis list is no longer necessary. */
	free (bslist);

	/* Sort the basis list according to parent boxes. */
	qsort (boxidx, nbs, 3 * sizeof(int), idxcomp);

	/* Count the number of non-empty boxes. */
	for (nebox = i = 1, idx = boxidx + 3; i < nbs; ++i, idx += 3) {
		/* Skip the current box if it was already counted. */
		if (!idxcomp (idx, idx - 3)) continue;

		++nebox;
	}

	/* Allocate the list of boxes. */
	boxlist = calloc (nebox, sizeof(boxdesc));

	/* Allocate the backend array. */
	rhsptr = rhsbuf = fftwf_malloc (nebox * nfftprod * sizeof(complex float));

	fprintf (stderr, "Rank %d: Expanded FFT buffer size size: %ld bytes\n",
			rank, nebox * nfftprod * sizeof(complex float));

	/* Set up the first box structure. */
	memcpy (boxlist[0].index, boxidx, 3 * sizeof(int));
	omp_init_lock (&(boxlist[0].lock));
	boxlist[0].rhs = rhsptr;
	rhsptr += nfftprod;

	/* Now loop through all bases and set up the box cache pointers. */
	for (nebox = i = 1, idx = boxidx + 3; i < nbs; ++i, idx += 3) {
		/* Skip a previously-counted box. */
		if (!idxcomp (idx, idx - 3)) continue;

		/* Set up the next box structure. */
		memcpy (boxlist[nebox].index, idx, 3 * sizeof(int));
		omp_init_lock (&(boxlist[nebox].lock));
		boxlist[nebox].rhs = rhsptr;
		rhsptr += nfftprod;

		/* Count this box. */
		++nebox;
	}

	/* Sort the box descriptors. */
	qsort (boxlist, nebox, sizeof(boxdesc), boxcomp);

	/* The box list is no longer necessary. */
	free (boxidx);

	return nebox;
}

/* Clear the cache arrays for a new iteration. */
void clrdircache () {
	int i;

	/* Blank the box fill marker. */
	for (i = 0; i < nebox; ++i) boxlist[i].fill = 0;
}

/* Free the allocated memory in the direct-interaction cache structure. */
void freedircache () {
	free (rhsbuf);
	free (boxlist);
	free (gridints);
}

/* Check the cache for the given RHS. If it exists, return the pre-cached copy.
 * Otherwise, cache the provided copy and take the DFT before returning a pointer
 * to the cache bin. */
complex float *cacheboxrhs (int *bslist, int nbs, int boxkey) {
	int l, i;
	boxdesc key, *lbox;
	complex float *bptr, *rhs;

	/* Get the index for the first basis in the box. */
	bsindex (bslist[0], key.index);

	/* The box index. */
	key.index[0] /= fmaconf.bspbox;
	key.index[1] /= fmaconf.bspbox;
	key.index[2] /= fmaconf.bspbox;

	/* Search for the box index. */
	lbox = bsearch (&key, boxlist, nebox, sizeof(boxdesc), boxcomp);

	/* The search failed for some reason. */
	if (!lbox) return NULL;

	/* Point to the RHS storage. */
	bptr = lbox->rhs;

	/* The cache check and fill operation must be thread safe. */
	omp_set_lock (&(lbox->lock));

	/* Cache miss. Fill the box. */
	if (!(lbox->fill)) {
		/* Clear the cache storage. */
		memset (bptr, 0, nfftprod * sizeof(complex float));

		/* Grab the local input vector for caching. */
		rhs = (complex float *)ScaleME_getInputVec (boxkey);
		
		/* Populate the local grid. */
		for (i = 0; i < nbs; ++i) {
			/* Find the basis index. */
			bsindex (bslist[i], key.index);
			
			/* Find the position in the local box. */
			key.index[0] %= fmaconf.bspbox;
			key.index[1] %= fmaconf.bspbox;
			key.index[2] %= fmaconf.bspbox;

			/* The index into the RHS array. */
			l = SQIDX(nfft[0],key.index[0],key.index[1],key.index[2]);
			
			/* Fill the expanded FFT grid. */
			bptr[l] = rhs[i];
		}
		
		/* Transform the cached RHS. */
		fftwf_execute_dft (fplan, bptr, bptr);
		
		/* Mark the cache spot as full. */
		lbox->fill = 1;

	}

	/* Free the lock to allow other threads to access the cache block. */
	omp_unset_lock (&(lbox->lock));

	/* Return the RHS. */
	return bptr;
}

/* Precompute some values for the direct interactions. */
int dirprecalc () {
	int totbpnbr, nborsvol, rank;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

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

	fprintf (stderr, "Rank %d: Green's function grid size: %ld bytes\n",
			rank, totbpnbr * sizeof(complex float));

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
void blockinteract (int tkey, int tct, int *skeys, int *scts, int numsrc) {
	int i, l, boxoff[3], idx[3], *obslist, *srclist;
	complex float *buf, *gptr, *bptr;
	complex float *cobs;

	/* Clear the local output buffer. */
	buf = fftwf_malloc (nfftprod * sizeof(complex float));
	memset (buf, 0, nfftprod * sizeof(complex float));

	/* Find the output vector segment and the target basis list. */
	obslist = ScaleME_getBasisList (tkey);
	cobs = (complex float *)ScaleME_getOutputVec (tkey);

	/* Find the index for the first basis in the target box. */
	bsindex (obslist[0], boxoff);

	/* Find the minimum box index for near interactions. */
	boxoff[0] = (boxoff[0] / fmaconf.bspbox) - fmaconf.numbuffer;
	boxoff[1] = (boxoff[1] / fmaconf.bspbox) - fmaconf.numbuffer;
	boxoff[2] = (boxoff[2] / fmaconf.bspbox) - fmaconf.numbuffer;

	/* Populate the local grid with contributions from near boxes. */
	for (l = 0; l < numsrc; ++l) {
		srclist = ScaleME_getBasisList (skeys[l]);

		/* Get the cached RHS for the source box in question.
		 * The cache may need to be filled. */
		bptr = cacheboxrhs (srclist, scts[l], skeys[l]);

		/* Get the index fo the first basis in the source box. */
		bsindex (srclist[0], idx);

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
	for (l = 0; l < tct; ++l) {
		bsindex (obslist[l], idx);

		/* Convert the global grid position to a local position. */
		idx[0] %= fmaconf.bspbox;
		idx[1] %= fmaconf.bspbox;
		idx[2] %= fmaconf.bspbox;

		/* Augment the RHS. */
		cobs[l] += buf[SQIDX(nfft[0],idx[0],idx[1],idx[2])];
	}

	free (buf);

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
				/* The self-term needs no special attention because
				 * the integration grids don't coincide. */
				*(grf++) = scale * rcvint (k0, zero, dist, cell, fsgreen);
			}
		}
	}

	return mex;
}

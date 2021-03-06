#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <math.h>
#include <fftw3.h>

/* Define OpenMP locking functions if locking is not used. */
#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_lock_t;
int omp_set_lock (omp_lock_t *x) { return 0; }
int omp_init_lock (omp_lock_t *x) { return 0; }
int omp_unset_lock (omp_lock_t *x) { return 0; }
#endif /* _OPENMP */

#include "ScaleME.h"

#include "precision.h"

#include "direct.h"
#include "util.h"
#include "fsgreen.h"
#include "mlfma.h"
#include "integrate.h"

/* Define macros for packed indexing. */
#define PPYR(i) (((i) * ((i) + 1) * ((i) + 2)) / 6)
#define PTRI(i) (((i) * ((i) + 1)) / 2)

typedef struct {
	int index[3], fill;
	omp_lock_t lock;
	cplx *rhs;
} boxdesc;

/* Buffers for the RHS cache, the Green's functions, and a workspace. */
static boxdesc *boxlist;
static cplx *gridints, *rhsbuf;
static int nbors, nfftprod, nfft, nebox;
static FFTW_PLAN fplan, bplan;

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
	cplx *rhsptr;
	int nbs, *bslist, i, *idx, *boxidx, rank;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* Get the complete list of locally required basis functions. */
	ScaleME_getLocallyReqBasis (&nbs, &bslist);

	/* Allocate the box index array. */
	boxidx = malloc (3 * nbs * sizeof(int));

	/* Build the array of box indices. */
	for (i = 0, idx = boxidx; i < nbs; ++i, idx += 3)
		GRID (idx, bslist[i], fmaconf.nx, fmaconf.ny);

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
	rhsptr = rhsbuf = FFTW_MALLOC (nebox * nfftprod * sizeof(cplx));

	fprintf (stderr, "Rank %d: Expanded FFT buffer size size: %ld bytes\n",
			rank, nebox * nfftprod * sizeof(cplx));

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
	FFTW_FREE (rhsbuf);
	FFTW_FREE (gridints);
	free (boxlist);
}

/* Check the cache for the given RHS. If it exists, return the pre-cached copy.
 * Otherwise, cache the provided copy and take the DFT before returning a pointer
 * to the cache bin. */
cplx *cacheboxrhs (int bsl, int boxkey) {
	int l, i;
	boxdesc key, *lbox;
	cplx *bptr, *rhs;

	/* Get the index for the first basis in the box. */
	GRID (key.index, bsl, fmaconf.nx, fmaconf.ny);

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
		memset (bptr, 0, nfftprod * sizeof(cplx));

		/* Grab the local input vector for caching. */
		rhs = (cplx *)ScaleME_getInputVec (boxkey);

		/* Populate the local grid. */
		for (i = 0; i < fmaconf.bspboxvol; ++i) {
			/* Find the position in the local box. */
			GRID(key.index, i, fmaconf.bspbox, fmaconf.bspbox);

			/* The index into the RHS array. */
			l = IDX(nfft,key.index[0],key.index[1],key.index[2]);

			/* Fill the expanded FFT grid. */
			bptr[l] = rhs[i];
		}

		/* Transform the cached RHS. */
		FFTW_EXECUTE_DFT (fplan, bptr, bptr);

		/* Mark the cache spot as full. */
		lbox->fill = 1;
	}

	/* Free the lock to allow other threads to access the cache block. */
	omp_unset_lock (&(lbox->lock));

	/* Return the RHS. */
	return bptr;
}

/* Precompute a cache of unique, integrated Green's function values in a packed
 * array to be later used to populate the near-field interactions. */
int greencache (cplx *grf, int n, real k0, real cell) {
	real zero[3] = { 0., 0., 0. }, dc[3] = { cell, cell, cell };

	/* Compute the self term. */
	grf[0] = rcvint (k0, NULL, zero, dc, duffyint, fsgrnduffy);

#pragma omp parallel default(shared)
{
	int i, j, k, m, l;
	real dist[3];

#pragma omp for
	for (l = 1; l < n; ++l) {
		/* Compute the first grid index from the linear index. */
		i = cbrt (6 * l);
		if (PPYR(i) > l) --i;

		/* Compute the remainder of the index. */
		m = l - PPYR(i);

		/* Compute the second grid index from the linear index remainder. */
		j = sqrt (2 * m);
		if (PTRI(j) > m) --j;

		/* Compute the final grid index. */
		k = m - PTRI(j);

		/* Compute the distance. */
		dist[0] = cell * i;
		dist[1] = cell * j;
		dist[2] = cell * k;

		/* Integrate the appropraite term. */
		grf[l] = rcvint (k0, zero, dist, dc, srcint, fsgreen);
	}
}

	return n;
}

/* Precompute some values for the direct interactions. */
int dirprecalc (int numsrcpts) {
	int totbpnbr, nborsvol, rank, sepmax, ncache;
	cplx *grc;

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* The number of neighbor boxes per dimension. */
	nbors = nborsvol = 2 * fmaconf.numbuffer + 1;
	/* The number of near-neighbor boxes total. */
	nborsvol *= nborsvol * nborsvol;

	/* The FFT size. */
	nfft = 2 * fmaconf.bspbox;
	nfftprod = nfft * nfft * nfft;

	/* Build the expanded grid. */
	totbpnbr = nfftprod * nborsvol;
	gridints = FFTW_MALLOC (totbpnbr * sizeof(cplx));

	fprintf (stderr, "Rank %d: Green's function grid size: %ld bytes\n",
			rank, totbpnbr * sizeof(cplx));

	/* The forward FFT plan transforms all boxes in one pass. */
	fplan = FFTW_PLAN_DFT_3D (nfft, nfft, nfft,
			gridints, gridints, FFTW_FORWARD, FFTW_MEASURE);
	/* The inverse FFT plan only transforms a single box. */
	bplan = FFTW_PLAN_DFT_3D (nfft, nfft, nfft,
			gridints, gridints, FFTW_BACKWARD, FFTW_MEASURE);

	/* This is the maximum single-index basis separation to be cached. */
	sepmax = nfft + fmaconf.numbuffer * fmaconf.bspbox;
	/* This is the size of the packed-storage value cache. */
	ncache = PPYR(sepmax);
	/* Allocate and fill the Green's function integration cache. */
	grc = malloc(ncache * sizeof(cplx));
	greencache (grc, ncache, fmaconf.k0, fmaconf.cell);
	if (!rank)
		fprintf (stderr, "Cached %d Green's function integrations\n", ncache);

#pragma omp parallel default(shared)
{
	int off[3], l, idx[3];
	cplx *grf;

#pragma omp for
	for (l = 0; l < nborsvol; ++l) {
		GRID(idx, l, nbors, nbors);

		grf = gridints + l * nfftprod;

		off[0] = (idx[0] - fmaconf.numbuffer) * fmaconf.bspbox;
		off[1] = (idx[1] - fmaconf.numbuffer) * fmaconf.bspbox;
		off[2] = (idx[2] - fmaconf.numbuffer) * fmaconf.bspbox;

		/* Build the Green's function grid for this local box. */
		greengrid (grf, fmaconf.bspbox, nfft, off, fmaconf.k0, grc);

		/* Fourier transform the Green's function. */
		FFTW_EXECUTE_DFT (fplan, grf, grf);
	}
}

	/* Allocate the local cache structure. */
	mkdircache ();

	/* Free the integration cache. */
	free (grc);

	return nfftprod;
}

/* Evaluate at a group of observers the fields due to a group of sources. */
void blockinteract (int tkey, int tct, int *skeys, int *scts, int numsrc) {
	int i, l, boxoff[3], idx[3], *obslist, *srclist;
	cplx *buf, *gptr, *bptr;
	cplx *cobs;

	/* Clear the local output buffer. */
	buf = FFTW_MALLOC (nfftprod * sizeof(cplx));
	memset (buf, 0, nfftprod * sizeof(cplx));

	/* Find the output vector segment and the target basis list. */
	obslist = ScaleME_getBasisList (tkey);
	cobs = (cplx *)ScaleME_getOutputVec (tkey);

	/* Find the index for the first basis in the target box. */
	GRID (boxoff, obslist[0], fmaconf.nx, fmaconf.ny);

	/* Find the minimum box index for near interactions. */
	boxoff[0] -= fmaconf.numbuffer;
	boxoff[1] -= fmaconf.numbuffer;
	boxoff[2] -= fmaconf.numbuffer;

	/* Populate the local grid with contributions from near boxes. */
	for (l = 0; l < numsrc; ++l) {
		srclist = ScaleME_getBasisList (skeys[l]);

		/* Get the cached RHS for the source box in question.
		 * The cache may need to be filled. */
		bptr = cacheboxrhs (srclist[0], skeys[l]);

		/* Get the index fo the first basis in the source box. */
		GRID (idx, srclist[0], fmaconf.nx, fmaconf.ny);

		/* Find the local box number. */
		idx[0] -= boxoff[0];
		idx[1] -= boxoff[1];
		idx[2] -= boxoff[2];

		/* Point to the Green's function for this box. */
		gptr = gridints + nfftprod * IDX(nbors,idx[0],idx[1],idx[2]);

		/* Convolve the source field with the Green's function and
		 * augment the field at the target. */
		for (i = 0; i < nfftprod; ++i) buf[i] += gptr[i] * bptr[i];
	}

	/* Inverse transform the grid in place. */
	FFTW_EXECUTE_DFT (bplan, buf, buf);

	/* Augment with output with the local convolution. */
	/* Note that each ScaleME "basis" is actually a finest-level group. */
	for (l = 0; l < fmaconf.bspboxvol; ++l) {
		GRID(idx, l, fmaconf.bspbox, fmaconf.bspbox);

		/* Augment the RHS. */
		cobs[l] += buf[IDX(nfft,idx[0],idx[1],idx[2])];
	}

	FFTW_FREE (buf);

	return;
}

/* Build the scaled, extended Green's function grf on an expanded cubic grid
 * with values taken from the cache grc. */
int greengrid (cplx *grf, int m, int mex, int *off, real k0, cplx *grc) {
	int ip, jp, kp, l, mt, idx[3];
	real scale;

	/* The total number of samples. */
	mt = mex * mex * mex;

	/* The scale of the integral equation solution. */
	scale = k0 * k0 / (real)mt;

	/* Compute the interactions. */
	for (l = 0; l < mt; ++l) {
		GRID(idx, l, mex, mex);

		/* Compute the source and observer separation in cell lengths. */
		ip = (idx[0] < m) ? idx[0] : (idx[0] - mex);
		ip = abs (ip - off[0]);

		jp = (idx[1] < m) ? idx[1] : (idx[1] - mex);
		jp = abs (jp - off[1]);

		kp = (idx[2] < m) ? idx[2] : (idx[2] - mex);
		kp = abs (kp - off[2]);

		/* The maximum separation index gets stored in idx[0]. */
		idx[0] = MAX(ip, MAX(jp, kp));
		/* The minimum separation index gets stored in idx[2]. */
		idx[2] = MIN(ip, MIN(jp, kp));
		/* The middle (remaining) index gets stored in idx[1]. */
		idx[1] = ip + jp + kp - idx[0] - idx[2];

		/* Copy the proper integration into place from the cache. */
		*(grf++) = scale * grc[PPYR(idx[0]) + PTRI(idx[1]) + idx[2]];
	}

	return mex;
}

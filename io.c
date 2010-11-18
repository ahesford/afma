#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "io.h"
#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"

#include "util.h"

/* Compare the basis map for sorting. */
int mapcomp (const void *left, const void *right) {
	long *il = (long *)left, *ir = (long *)right;
	return *il - *ir;
}

/* Skip comments in the input file. */
void skipcomments (FILE *fp) {
	fpos_t loc;
	char buf[1024];

	do {
		fgetpos (fp, &loc);
		if (!fgets (buf, 1024, fp)) break;
	} while (strstr (buf, "#") != NULL);

	fsetpos (fp, &loc);
}

/* Read the DBIM configuration. */
void getdbimcfg (char *fname, int *maxit, float *regparm, float *tol) {
	FILE *fp;
	char buf[1024];

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return;
	}

	/* Read the number of DBIM iterations. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (sscanf (buf, "%d %d", maxit, maxit + 1) < 2) maxit[1] = maxit[0];

	/* Read the regularization bounds, step and skip. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	/* The skip defaults to 1, if none is provided. */
	if (sscanf (buf, "%f %f %f %f", regparm, regparm + 1,
				regparm + 2, regparm + 3) < 4)
		regparm[3] = 1;

	/* Read the DBIM tolerance. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (sscanf (buf, "%f %f", tol, tol + 1) < 2) tol[1] = tol[2] / 10;

	fclose (fp);
}

/* Read the configuration file and set parameters. */
void getconfig (char *fname, solveparm *hislv, solveparm *loslv,
		measdesc *src, measdesc *obs, int obscount) {
	FILE *fp;
	char buf[1024];
	int nmax, nbox, i;

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return;
	}

	/* Read the number of cells. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d %d", &(fmaconf.nx), &(fmaconf.ny), &(fmaconf.nz));

	/* Read the lower corner of the domain. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %f %f", fmaconf.cen, fmaconf.cen + 1,
			fmaconf.cen + 2, &(fmaconf.cell));

	fmaconf.min[0] = fmaconf.cen[0] - 0.5 * (float)(fmaconf.nx) * fmaconf.cell;
	fmaconf.min[1] = fmaconf.cen[1] - 0.5 * (float)(fmaconf.ny) * fmaconf.cell;
	fmaconf.min[2] = fmaconf.cen[2] - 0.5 * (float)(fmaconf.nz) * fmaconf.cell;

	fmaconf.cellvol = fmaconf.cell * fmaconf.cell * fmaconf.cell;

	/* Set the wave number to 2 pi, since wavelength is the length unit. */
	fmaconf.k0 = 2 * M_PI;

	/* Read the MLFMA level and fast translation configuration. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (sscanf (buf, "%d %d %d %d %d", &(fmaconf.bspbox), 
				&(fmaconf.toplev), &(fmaconf.fo2iterm),
				&(fmaconf.fo2iord), &(fmaconf.fo2iosr)) < 5)
		fmaconf.fo2iterm = fmaconf.fo2iord = fmaconf.fo2iosr = 0;

	/* The total number of basis functions in a box. */
	fmaconf.bspboxvol = fmaconf.bspbox * fmaconf.bspbox * fmaconf.bspbox;

	/* Compute the number of finest-level groups for use with the FMM. */
	fmaconf.nx = GEDIV(fmaconf.nx,fmaconf.bspbox);
	fmaconf.ny = GEDIV(fmaconf.ny,fmaconf.bspbox);
	fmaconf.nz = GEDIV(fmaconf.nz,fmaconf.bspbox);

	/* The maximum FMM length. */
	nmax = MAX(fmaconf.nx,MAX(fmaconf.ny,fmaconf.nz));

	/* Set the global number of FMM bases, for easy reference. */
	fmaconf.gnumbases = fmaconf.nx * fmaconf.ny * fmaconf.nz;

	/* The length of a finest-level group. */
	fmaconf.grplen = fmaconf.bspbox * fmaconf.cell;

	/* Compute the maximum FMM level for the desired finest level density. */
	for (fmaconf.maxlev = 2, nbox = 4; nbox < nmax; nbox <<= 1, ++fmaconf.maxlev);

	/* Read the number of MLFMA buffer boxes. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.numbuffer));

	/* Read the MLFMA precision. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(fmaconf.precision));

	/* Read the MLFMA interpolation order. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.interpord));

	/* Read the high-accuracy iterative solver configuration. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d %f", &(hislv->maxit),
			&(hislv->restart), &(hislv->epscg));

	/* Read the low-accuracy iterative solver configuration, if desired. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (loslv)
		sscanf (buf, "%d %d %f", &(loslv->maxit),
				&(loslv->restart), &(loslv->epscg));

	/* Read the source radius and ignore it since plane waves are used. */
	skipcomments (fp);
	fgets (buf, 1024, fp);

	/* Read the source theta values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", src->trange, src->trange + 1, &(src->ntheta));

	/* Convert the degree values to radians. */
	src->trange[0] *= M_PI / 180;
	src->trange[1] *= M_PI / 180;

	/* Read the source phi values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", src->prange, src->prange + 1, &(src->nphi));

	/* Convert the degree values to radians. */
	src->prange[0] *= M_PI / 180;
	src->prange[1] *= M_PI / 180;

	/* Read the observer radius and skip it since plane waves are used. */
	skipcomments (fp);
	fgets (buf, 1024, fp);

	/* Read the specified number of observation configurations. */
	for (i = 0; i < obscount && !feof(fp); ++i) {
		/* Read the observer theta values. */
		skipcomments (fp);
		fgets (buf, 1024, fp);
		sscanf (buf, "%f %f %d", obs[i].trange, obs[i].trange + 1, &(obs[i].ntheta));
		
		/* Convert the degree values to radians. */
		obs[i].trange[0] *= M_PI / 180;
		obs[i].trange[1] *= M_PI / 180;
		
		/* Read the observer phi values. */
		skipcomments (fp);
		fgets (buf, 1024, fp);
		sscanf (buf, "%f %f %d", obs[i].prange, obs[i].prange + 1, &(obs[i].nphi));
		
		/* Convert the degree values to radians. */
		obs[i].prange[0] *= M_PI / 180;
		obs[i].prange[1] *= M_PI / 180;
	}

	fclose (fp);
}

/* Read the gridded file into local arrays. */
int getcontrast (complex float *contrast, char *fname, int *size,
		int *bslist, int nbs, int bpb) {
	long i, *bsmap, l, nelt, nnz;
	int k, bsize[3], idx[3], gidx[3], bpbvol = bpb * bpb * bpb;

	/* MPI data types for file I/O. */
	MPI_Status stat;
	MPI_Datatype cplx;
	MPI_File fh;

	nelt = (long)nbs * (long)bpbvol;

	/* Zero out the contrast in case nothing can be read. */
	memset (contrast, 0, nelt * sizeof(complex float));

	/* Define and commit a complex type for MPI reads. */
	MPI_Type_contiguous (2, MPI_FLOAT, &cplx);
	MPI_Type_commit (&cplx);

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_RDONLY,
				MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return 0;
	}

	/* Read the basis grid size. */
	MPI_File_read (fh, bsize, 3, MPI_INT, &stat);

	/* Position the file to read complex values after the header. */
	MPI_File_set_view (fh, 3 * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);

	/* Allocate the map from the read order to the local storage order. */
	bsmap = malloc (2L * nelt * sizeof(long));

	/* Build a basis map for sorting. */
	for (nnz = i = l = 0; i < nbs; ++i) {
		bsindex (bslist[i], gidx);

		/* Find the starting basis index from the group index. */
		gidx[0] *= bpb;
		gidx[1] *= bpb;
		gidx[2] *= bpb;

		for (k = 0; k < bpbvol; ++k, ++l) {
			/* Compute the index of the basis within the group. */
			GRID(idx, bpb, k);
			idx[0] += gidx[0];
			idx[1] += gidx[1];
			idx[2] += gidx[2];

			/* Don't read this element if it exceeds grid bounds. */
			if (idx[0] >= bsize[0] || idx[1] >= bsize[1] || 
					idx[2] >= bsize[2]) continue;

			/* Compute the basis position. */
			bsmap[2 * nnz] = (long)idx[0] +
				(long)bsize[0] * (long)idx[1] + 
				(long)bsize[0] * (long)bsize[1] * (long)idx[2];
			bsmap[2 * nnz + 1] = l;
			++nnz;
		}
	}

	/* Sort the nonzero basis map entries. */
	qsort (bsmap, nnz, 2 * sizeof(long), mapcomp);

	/* Read the nonzero elements. */
	for (i = 0; i < nnz; ++i)
		MPI_File_read_at (fh, bsmap[2*i], contrast+bsmap[2*i+1], 1, cplx, &stat);

	/* Free the sorted map and MPI complex type, and close the file. */
	free (bsmap);
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	return nbs;
}

/* Distributed write of a gridded file into local arrays. */
int prtcontrast (char *fname, complex float *crt, int *size,
		int *bslist, int nbs, int bpb) {
	int k, mpirank, bpbvol = bpb * bpb * bpb, bsize[3], gidx[3], idx[3];
	long i, *bsmap, nelt, l;
	
	/* MPI data types for file I/O. */
	MPI_Offset offset;
	MPI_Status stat;
	MPI_Datatype cplx;
	MPI_File fh;

	nelt = (long)nbs * (long)bpbvol;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Compute the output file size and store the size array. */
	offset = (long)size[0] * (long)size[1] * (long)size[2] * 
		(long)bpbvol * sizeof(complex float) + 3 * sizeof(int);

	bsize[0] = bpb * size[0];
	bsize[1] = bpb * size[1];
	bsize[2] = bpb * size[2];

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
				MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: could not open %s.\n", fname);
		return 0;
	}

	/* Set the size of the output file and set the view to write the header. */
	MPI_File_set_size (fh, offset);
	MPI_File_set_view (fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	/* The first process should write the size header in the file. */
	if (!mpirank) MPI_File_write (fh, bsize, 3, MPI_INT, &stat);
	/* Ensure the write is committed and all processes should wait. */
	MPI_File_sync (fh);
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_File_sync (fh);

	MPI_Type_contiguous (2, MPI_FLOAT, &cplx);
	MPI_Type_commit (&cplx);

	/* Set the file view to skip over the grid size header. */
	MPI_File_set_view (fh, 3 * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);

	/* Allocate the map from the read order to the local storage order. */
	bsmap = malloc (2 * nelt * sizeof(long));

	/* Build a basis map for sorting. */
	for (i = l = 0; i < nbs; ++i) {
		bsindex (bslist[i], gidx);

		/* Find the starting basisindex from the group index. */
		gidx[0] *= bpb;
		gidx[1] *= bpb;
		gidx[2] *= bpb;

		for (k = 0; k < bpbvol; ++k, ++l) {
			/* Compute the index of the basis within the group. */
			GRID(idx, bpb, k);
			idx[0] += gidx[0];
			idx[1] += gidx[1];
			idx[2] += gidx[2];

			/* Compute the basis position. */
			bsmap[2 * l] = (long)idx[0] +
				(long)bsize[0] * (long)idx[1] + 
				(long)bsize[0] * (long)bsize[1] * (long)idx[2];
			bsmap[2 * l + 1] = l;
		}
	}

	/* Sort the basis map. */
	qsort (bsmap, nelt, 2 * sizeof(long), mapcomp);

	/* Perform an ordered write of the file. */
	for (i = 0; i < nelt; ++i)
		MPI_File_write_at(fh, bsmap[2*i], crt+bsmap[2*i+1], 1, cplx, &stat);

	/* Free the basis map. */
	free (bsmap);

	/* Free the MPI file type and close the file. */
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	return size[0] * size[1] * size[2];
}

/* Append the provided field to the specified field file. */
int writefld (char *fname, int nrows, int ncols, complex float *field) {
	FILE *fp;
	int size[2] = { nrows, ncols }, count = nrows * ncols;

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not write file %s.\n", fname);
		return 0;
	}

	/* Write the matrix size. */
	fwrite (size, sizeof(int), 2, fp);
	/* Write the values. */
	fwrite (field, sizeof(complex float), count, fp);
	fclose (fp);

	return count;
}

int readfld (complex float *field, char *fname, int nobs) {
	FILE *fp;
	int size[2];

	if (!(fp = fopen(fname, "r"))) {
		fprintf (stderr, "ERROR: could not read file %s.\n", fname);
		return 0;
	}

	/* Read the stored matrix size. */
	fread (size, sizeof(int), 2, fp);

	if (size[0] * size[1] != nobs) {
		fprintf (stderr, "ERROR: wrong matrix size in file %s.\n", fname);
		fclose (fp);
		return 0;
	}

	/* Read the values. */
	fread (field, sizeof(complex float), nobs, fp);
	fclose (fp);

	return nobs;
}

int getfields (char *inproj, complex float *field, int nobs, int nsrc, float *nrm) {
	complex float *fldptr;
	int i, j;
	float lerr;
	char fname[1024], fmt[1024];

	/* Find the width of the integer label in the field name. */
	j = (int)ceil(log10(nsrc));
	sprintf (fmt, "%%s.tx%%0%dd.field", j);

	/* The initial norm of the matrix is zero, if it is desired. */
	if (nrm) *nrm = 0;

	/* Read each column of the matrix and compute the norm. */
	for (fldptr = field, i = 0; i < nsrc; ++i, fldptr += nobs) {
		sprintf (fname, fmt, inproj, i);
		readfld (fldptr, fname, nobs);

		/* If the norm isn't desired, don't comput anything else. */
		if (!nrm) continue;
		
		for (j = 0; j < nobs; ++j) {
			lerr = cabs (fldptr[j]);
			*nrm += lerr * lerr;
		}
	}

	return nobs * nsrc;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "precision.h"

#include "io.h"
#include "util.h"

#define HLEN 5

typedef struct {
	int bsi;
	long off;
} ctmap;

/* Compare the map entries for sorting. */
static int ctmapcomp (const void *left, const void *right) {
	return ((ctmap *)left)->bsi - ((ctmap *)right)->bsi;
}

/* Sort the contrast file according to the basis map, which will be constructed. */
static int sortbsmap (ctmap *map, int *bsl, int nbs, int bpvol) {
	int i;
	long l;

	/* Build the contrast map. */
	for (l = i = 0; i < nbs; ++i, l += bpvol) {
		map[i].bsi = bsl[i];
		map[i].off = l;
	}

	/* Sort the contrast map. */
	qsort (map, nbs, sizeof(ctmap), ctmapcomp);

	return nbs;
}

/* Read the gridded file into local arrays. */
int getctgrp (cplx *crt, char *fname, int *size, int *bsl, int nbs, int bpb) {
	long i, nelt;
	int bsize[HLEN], bpbvol = bpb * bpb * bpb;
	ctmap *map;

	MPI_File fh;
	MPI_Datatype cplx;
	MPI_Status stat;

	nelt = (long)nbs * (long)bpbvol;

	/* Zero out the contrast in case nothing can be read. */
	memset (crt, 0, nelt * sizeof(cplx));

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_RDONLY,
				MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return 0;
	}

	/* Read the basis grid size. */
	MPI_File_set_view (fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_read (fh, bsize, HLEN, MPI_INT, &stat);

	if (bsize[0] != 0) {
		fprintf (stderr, "ERROR: %s is not a group-ordered file.\n", fname);
		MPI_File_close (&fh);
		return 0;
	} else if (bsize[1] != size[0] || bsize[2] != size[1] || bsize[3] != size[2]) {
		fprintf (stderr, "ERROR: %s does not match group count.\n", fname);
		MPI_File_close (&fh);
		return 0;
	} else if (bsize[4] != bpb) {
		fprintf (stderr, "ERROR: %s does not match group size.\n", fname);
		MPI_File_close (&fh);
		return 0;
	}

	/* Prepare the basis sorting map. */
	map = malloc (nbs * sizeof(ctmap));
	sortbsmap (map, bsl, nbs, bpbvol);

	/* Create the datatype storing all values in a group. */
	MPI_Type_contiguous (2 * bpbvol, MPIREAL, &cplx);
	MPI_Type_commit (&cplx);

	/* Perform the ordered, grouped read from the file. */
	MPI_File_set_view (fh, HLEN * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);
	for (i = 0; i < nbs; ++i)
		MPI_File_read_at (fh, map[i].bsi, crt + map[i].off, 1, cplx, &stat);

	/* Fre the MPI file type and close the file. */
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	free (map);

	return nbs;
}

/* Distributed write of a gridded file into local arrays. */
int prtctgrp (char *fname, cplx *crt, int *size, int *bsl, int nbs, int bpb) {
	int mpirank, bpbvol = bpb * bpb * bpb, bhdr[HLEN];
	long i, nelt;
	ctmap *map;

	MPI_File fh;
	MPI_Datatype cplx;
	MPI_Status stat;

	nelt = (long)nbs * (long)bpbvol;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_WRONLY |
				MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: could not open %s.\n", fname);
		return 0;
	}

	/* Re-sort the contrast map according to the group indices. */
	map = malloc (nbs * sizeof(ctmap));
	sortbsmap (map, bsl, nbs, bpbvol);

	/* Create the datatype storing all values in a group. */
	MPI_Type_contiguous (2 * bpbvol, MPIREAL, &cplx);
	MPI_Type_commit (&cplx);

	/* Build the special header to note the file is group-ordered. */
	bhdr[0] = 0;
	bhdr[1] = size[0];
	bhdr[2] = size[1];
	bhdr[3] = size[2];
	bhdr[4] = bpb;

	/* Write the file header, if the process is the root. */
	MPI_File_set_view (fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	if (!mpirank) MPI_File_write (fh, bhdr, HLEN, MPI_INT, &stat);

	/* Write the values from the buffer in one pass. */
	MPI_File_set_view (fh, HLEN * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);
	for (i = 0; i < nbs; ++i)
		MPI_File_write_at (fh, map[i].bsi, crt + map[i].off, 1, cplx, &stat);

	/* Free the MPI file type and close the file. */
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	free (map);

	return nbs;
}

/* Append the provided field to the specified field file. */
int writefld (char *fname, int nrows, int ncols, cplx *field) {
	FILE *fp;
	int size[2] = { nrows, ncols }, count = nrows * ncols;

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not write file %s.\n", fname);
		return 0;
	}

	/* Write the matrix size. */
	fwrite (size, sizeof(int), 2, fp);
	/* Write the values. */
	fwrite (field, sizeof(cplx), count, fp);
	fclose (fp);

	return count;
}

int readfld (cplx *field, char *fname, int nobs) {
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
	fread (field, sizeof(cplx), nobs, fp);
	fclose (fp);

	return nobs;
}

int getfields (char *inproj, cplx *field, int nobs, int nsrc, real *nrm) {
	cplx *fldptr;
	int i, j;
	real lerr;
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

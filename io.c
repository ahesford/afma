#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "io.h"
#include "util.h"

/* Read the gridded file into local arrays. */
int getcontrast (complex float *contrast, char *fname, int *size,
		int *bslist, int nbs, int bpb) {
	long i, l, nelt, offset;
	int k, bsize[3], idx[3], gidx[3], bpbvol = bpb * bpb * bpb, rrow;
	FILE *fh;

	nelt = (long)nbs * (long)bpbvol;

	/* Zero out the contrast in case nothing can be read. */
	memset (contrast, 0, nelt * sizeof(complex float));

	/* Open the MPI file with the specified name. */
	if (!(fh = fopen(fname, "rb"))) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return 0;
	}

	/* Read the basis grid size. */
	fread (bsize, sizeof(int), 3, fh);

	/* Build a basis map for sorting. */
	for (i = l = 0; i < nbs; ++i) {
		GRID (gidx, bslist[i], size[0], size[1]);

		/* Find the starting basis index from the group index. */
		gidx[0] *= bpb;
		gidx[1] *= bpb;
		gidx[2] *= bpb;

		/* Don't try to read beyond the bounds of the file. */
		rrow = MIN(bpb, MAX(0, bsize[0] - gidx[0]));

		for (k = 0; k < bpbvol; k += bpb, l += bpb) {
			/* Compute the index of the basis within the group. */
			GRID(idx, k, bpb, bpb);
			idx[0] += gidx[0];
			idx[1] += gidx[1];
			idx[2] += gidx[2];

			/* Don't read this row if it exceeds grid bounds. */
			if (idx[1] >= bsize[1] || idx[2] >= bsize[2]) continue;

			/* Compute the basis position. */
			offset = 3L * sizeof(int) + sizeof(complex float) *
				((long)idx[0] + (long)bsize[0] * (long)idx[1] + 
				 (long)bsize[0] * (long)bsize[1] * (long)idx[2]);
			fseek (fh, offset, SEEK_SET);
			fread (contrast + l, sizeof(complex float), rrow, fh);
		}
	}

	fclose (fh);
	return nbs;
}

/* Distributed write of a gridded file into local arrays. */
int prtcontrast (char *fname, complex float *crt, int *size,
		int *bslist, int nbs, int bpb) {
	int fd, k, mpirank, bpbvol = bpb * bpb * bpb, bsize[3], gidx[3], idx[3];
	long i, nelt, l, offset;
	FILE *fh;
	
	nelt = (long)nbs * (long)bpbvol;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Compute the output file size and store the size array. */
	offset = (long)size[0] * (long)size[1] * (long)size[2] * 
		(long)bpbvol * sizeof(complex float) + 3 * sizeof(int);

	bsize[0] = bpb * size[0];
	bsize[1] = bpb * size[1];
	bsize[2] = bpb * size[2];

	/* Open the MPI file with the specified name. */
	if (!(fh = fopen(fname, "wb"))) {
		fprintf (stderr, "ERROR: could not open %s.\n", fname);
		return 0;
	}

	/* Set the size of the output file. */
	fd = fileno (fh);
	ftruncate (fd, offset);
	fsync (fd);
	MPI_Barrier (MPI_COMM_WORLD);
	fsync (fd);

	/* The first process should write the size header in the file. */
	if (!mpirank) fwrite (bsize, sizeof(int), 3, fh);

	/* Build a basis map for sorting. */
	for (i = l = 0; i < nbs; ++i) {
		GRID (gidx, bslist[i], size[0], size[1]);

		/* Find the starting basisindex from the group index. */
		gidx[0] *= bpb;
		gidx[1] *= bpb;
		gidx[2] *= bpb;

		for (k = 0; k < bpbvol; k += bpb, l += bpb) {
			/* Compute the index of the basis within the group. */
			GRID(idx, k, bpb, bpb);
			idx[0] += gidx[0];
			idx[1] += gidx[1];
			idx[2] += gidx[2];

			/* Compute the basis position. */
			offset = 3L * sizeof(int) + sizeof(complex float) * 
				((long)idx[0] + (long)bsize[0] * (long)idx[1] + 
				(long)bsize[0] * (long)bsize[1] * (long)idx[2]);
			fseek (fh, offset, SEEK_SET);
			fwrite (crt + l, sizeof(complex float), bpb,  fh);
		}
	}

	/* Free the MPI file type and close the file. */
	fclose (fh);
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

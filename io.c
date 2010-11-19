#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "io.h"
#include "util.h"

/* Read the gridded file into local arrays. */
int getctgrp (complex float *contrast, char *fname, int *size,
		int *bslist, int nbs, int bpb) {
	long i, l, nelt, offset;
	int k, bsize[5], idx[3], gidx[3], bpbvol = bpb * bpb * bpb, rrow;
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
	fread (bsize, sizeof(int), 5, fh);

	if (bsize[0] != 0) {
		fprintf (stderr, "ERROR: %s is not a group-ordered file.\n", fname);
		return 0;
	} else if (bsize[1] != size[0] || bsize[2] != size[1] || bsize[3] != size[2]) {
		fprintf (stderr, "ERROR: %s does not match group count.\n", fname);
		return 0;
	} else if (bsize[4] != bpb) {
		fprintf (stderr, "ERROR: %s does not match group size.\n", fname);
		return 0;
	}
	
	/* Perform the new, grouped read from the file. */
	for (l = i = 0; i < nbs; ++i, l += bpbvol) {
		/* Find the offset into the file. */
		offset = 5L * sizeof(int) + 
			(long)bslist[i] * (long)bpbvol * sizeof(complex float);
		fseek (fh, offset, SEEK_SET);
		fread (contrast + l, sizeof(complex float), bpbvol, fh);
	}

	fclose (fh);
	return nbs;
}

/* Distributed write of a gridded file into local arrays. */
int prtctgrp (char *fname, complex float *crt, int *size,
		int *bslist, int nbs, int bpb) {
	int mpirank, bpbvol = bpb * bpb * bpb, bhdr[5];
	long i, nelt, l;

	MPI_File fh;
	MPI_Datatype cplxgrp;
	MPI_Status stat;
	MPI_Offset offset;
	
	nelt = (long)nbs * (long)bpbvol;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Compute the output file size and store the size array. */
	offset = (long)size[0] * (long)size[1] * (long)size[2] * 
		(long)bpbvol * sizeof(complex float) + 3 * sizeof(int);

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | 
				MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: could not open %s.\n", fname);
		return 0;
	}

	/* Set the size of the output file. */
	MPI_File_set_size (fh, offset);
	MPI_File_set_view (fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

	/* Write the special header to note the file is group-ordered. */
	bhdr[0] = 0;
	bhdr[1] = size[0];
	bhdr[2] = size[1];
	bhdr[3] = size[2];
	bhdr[4] = bpb;

	/* The first process should write the size header in the file. */
	if (!mpirank) MPI_File_write (fh, bhdr, 5, MPI_INT, &stat);

	/* Ensure the write is committed. */
	MPI_File_sync (fh);
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_File_sync (fh);

	/* Create the datatype storing all values in a group. */
	MPI_Type_contiguous (2 * bpbvol, MPI_FLOAT, &cplxgrp);
	MPI_Type_commit (&cplxgrp);

	/* Set the file view to skip the header and point to groups. */
	MPI_File_set_view (fh, 5 * sizeof(int), cplxgrp, cplxgrp, "native", MPI_INFO_NULL);

	/* Write the values for each group in one pass. */
	for (i = l = 0; i < nbs; ++i, l += bpbvol) 
		MPI_File_write_at (fh, bslist[i], crt + l, 1, cplxgrp, &stat);

	/* Free the MPI file type and close the file. */
	MPI_File_close (&fh);
	MPI_Type_free (&cplxgrp);
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

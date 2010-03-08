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
		measdesc *src, measdesc *obs) {
	FILE *fp;
	char buf[1024];
	int nmax, nbox;

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return;
	}

	/* Read the number of cells. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d %d", &(fmaconf.nx), &(fmaconf.ny), &(fmaconf.nz));

	nmax = MAX(fmaconf.nx,MAX(fmaconf.ny,fmaconf.nz));

	/* Read the lower corner of the domain. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %f %f", fmaconf.cen, fmaconf.cen + 1,
			fmaconf.cen + 2, &(fmaconf.cell));

	fmaconf.min[0] = fmaconf.cen[0] - 0.5 * (float)(fmaconf.nx) * fmaconf.cell;
	fmaconf.min[1] = fmaconf.cen[1] - 0.5 * (float)(fmaconf.ny) * fmaconf.cell;
	fmaconf.min[2] = fmaconf.cen[2] - 0.5 * (float)(fmaconf.nz) * fmaconf.cell;

	fmaconf.cellvol = fmaconf.cell * fmaconf.cell * fmaconf.cell;

	/* Set the global number of bases, for easy reference. */
	fmaconf.gnumbases = fmaconf.nx * fmaconf.ny * fmaconf.nz;

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

	/* Compute the maximum FMM level for the desired finest level density. */
	for (fmaconf.maxlev = 0, nbox = fmaconf.bspbox; nbox < nmax; nbox <<= 1, ++fmaconf.maxlev);

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

	/* Read the source radius. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(src->radius));

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

	/* Read the observer radius. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(obs->radius));

	/* Read the observer theta values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", obs->trange, obs->trange + 1, &(obs->ntheta));

	/* Convert the degree values to radians. */
	obs->trange[0] *= M_PI / 180;
	obs->trange[1] *= M_PI / 180;

	/* Read the observer phi values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", obs->prange, obs->prange + 1, &(obs->nphi));

	/* Convert the degree values to radians. */
	obs->prange[0] *= M_PI / 180;
	obs->prange[1] *= M_PI / 180;

	fclose (fp);
}

/* Read the gridded file into local arrays. */
int getcontrast (complex float *contrast, char *fname, int *bslist, int nbs) {
	long i, *bsmap;

	/* MPI data types for file I/O. */
	MPI_Status stat;
	MPI_Datatype cplx;
	MPI_File fh;

	/* Zero out the contrast in case nothing can be read. */
	memset (contrast, 0, nbs * sizeof(complex float));

	/* Define and commit a complex type for MPI reads. */
	MPI_Type_contiguous (2, MPI_FLOAT, &cplx);
	MPI_Type_commit (&cplx);

	/* Open the MPI file with the specified name. */
	if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_RDONLY,
				MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return 0;
	}

	MPI_File_set_view (fh, 3 * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);

	/* Allocate the map from the read order to the local storage order. */
	bsmap = malloc (2 * nbs * sizeof(long));

	/* Build a basis map for sorting. */
	for (i = 0; i < nbs; ++i) {
		bsmap[2 * i] = bslist[i];
		bsmap[2 * i + 1] = i;
	}

	/* Sort the basis map. */
	qsort (bsmap, nbs, 2 * sizeof(long), mapcomp);

	for (i = 0; i < nbs; ++i)
		MPI_File_read_at (fh, bsmap[2*i], contrast+bsmap[2*i+1], 1, cplx, &stat);

	/* Free the sorted map and MPI complex type, and close the file. */
	free (bsmap);
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	return nbs;
}

/* Distributed write of a gridded file into local arrays. */
int prtcontrast (char *fname, complex float *crt, int *size, int *bslist, int nbs) {
	int mpirank;
	long i, *bsmap;
	
	/* MPI data types for file I/O. */
	MPI_Offset offset;
	MPI_Status stat;
	MPI_Datatype cplx;
	MPI_File fh;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Compute the output file size and store the size array. */
	offset = size[0] * size[1] * size[2] * sizeof(complex float) + 3 * sizeof(int);

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
	if (!mpirank) MPI_File_write (fh, size, 3, MPI_INT, &stat);
	/* Ensure the write is committed and all processes should wait. */
	MPI_File_sync (fh);
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_File_sync (fh);

	MPI_Type_contiguous (2, MPI_FLOAT, &cplx);
	MPI_Type_commit (&cplx);

	/* Set the file view to skip over the grid size header. */
	MPI_File_set_view (fh, 3 * sizeof(int), cplx, cplx, "native", MPI_INFO_NULL);

	/* Allocate the map from the read order to the local storage order. */
	bsmap = malloc (2 * nbs * sizeof(long));

	/* Build a basis map for sorting. */
	for (i = 0; i < nbs; ++i) {
		bsmap[2 * i] = bslist[i];
		bsmap[2 * i + 1] = i;
	}

	/* Sort the basis map. */
	qsort (bsmap, nbs, 2 * sizeof(long), mapcomp);

	/* Perform an ordered write of the file. */
	for (i = 0; i < nbs; ++i)
		MPI_File_write_at(fh, bsmap[2*i], crt+bsmap[2*i+1], 1, cplx, &stat);

	/* Free the basis map. */
	free (bsmap);

	/* Free the MPI file type and close the file. */
	MPI_File_close (&fh);
	MPI_Type_free (&cplx);

	return size[0] * size[1] * size[2];
}

/* Append the provided field to the specified field file. */
int writefld (char *fname, measdesc *obs, complex float *field) {
	FILE *fp;
	int size[2];

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not write file %s.\n", fname);
		return 0;
	}

	size[0] = obs->nphi;
	size[1] = obs->ntheta;

	/* Write the matrix size. */
	fwrite (size, sizeof(int), 2, fp);
	/* Write the values. */
	fwrite (field, sizeof(complex float), obs->count, fp);
	fclose (fp);

	return obs->count;
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
	char fname[1024];

	/* The initial norm of the matrix is zero, if it is desired. */
	if (nrm) *nrm = 0;

	/* Read each column of the matrix and compute the norm. */
	for (fldptr = field, i = 0; i < nsrc; ++i, fldptr += nobs) {
		sprintf (fname, "%s.%d.field", inproj, i);
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

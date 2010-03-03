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

void skipcomments (FILE *fp) {
	fpos_t loc;
	char buf[1024];

	do {
		fgetpos (fp, &loc);
		if (!fgets (buf, 1024, fp)) break;
	} while (strstr (buf, "#") != NULL);

	fsetpos (fp, &loc);
}

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

/* Read a portion of the contrast file and store it. */
void getcontrast (complex float *contrast, char *fname, int *bslist, int nbs) {
	int i, err;
	long offset;

	MPI_File fh;
	MPI_Offset spos;
	MPI_Status stat;

	/* Zero out the contrast in case nothing can be read. */
	memset (contrast, 0, nbs * sizeof(complex float));

	/* Open the file using MPI. */
	err = MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	if (err != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return;
	}

	/* Skip the contrast grid size. */
	MPI_File_seek (fh, 3 * sizeof(int), MPI_SEEK_SET);

	MPI_Barrier (MPI_COMM_WORLD);

	/* Read all values required by this processor. */
	for (offset = 0, i = 0; i < nbs; ++i) {
		spos = (long)(bslist[i] - offset) * sizeof(complex float);
		offset = bslist[i] + 1;

		/* Seek to the place holding the contrast. */
		MPI_File_seek (fh, spos, MPI_SEEK_CUR);
		/* Read the contrast value. */
		MPI_File_read (fh, contrast + i, 2, MPI_FLOAT, &stat);
	}

	MPI_File_close (&fh);
}

/* Write the contrast values using MPI files to keep things localized. */
int prtcontrast (char *fname, complex float *contrast,
		int *bslist, int nbs, int size[3]) {
	int i, err, mpirank;
	long offset;

	MPI_File fh;
	MPI_Offset spos;
	MPI_Status stat;

	/* Grab the process rank. */
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	/* Create and open the file using MPI. */
	err = MPI_File_open (MPI_COMM_WORLD, fname,
			MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

	if (err != MPI_SUCCESS) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return 0;
	}

	/* Write the grid size to the start of the file. */
	if (!mpirank) MPI_File_write (fh, size, 3, MPI_INT, &stat);

	/* Wait for the root node to write the size. */
	MPI_Barrier (MPI_COMM_WORLD);

	/* Skip over the grid size. */
	MPI_File_seek (fh, 3 * sizeof(int), MPI_SEEK_SET);

	/* Write all values on this processor. */
	for (offset = 0, i = 0; i < fmaconf.numbases; ++i) {
		spos = (long)(bslist[i] - offset) * sizeof(complex float);
		offset = bslist[i] + 1;

		/* Seek to the place to write the current contrast. */
		MPI_File_seek (fh, spos, MPI_SEEK_CUR);
		/* Write the contrast value. */
		MPI_File_write (fh, contrast + i, 2, MPI_FLOAT, &stat);
	}

	MPI_File_close (&fh);

	return fmaconf.gnumbases;
}

/* Print the header for the overall field matrix. */
int prtfldhdr (char *fname, measdesc *src, measdesc *obs) {
	FILE *fp;
	int size[2];

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not open field output.\n");
		return 0;
	}

	size[0] = obs->count;
	size[1] = src->count;

	fwrite (size, sizeof(int), 2, fp);

	fclose (fp);

	return size[0] * size[1];
}

/* Append the provided field to the specified field file. */
int appendfld (char *fname, measdesc *obs, complex float *field) {
	FILE *fp;

	if (!(fp = fopen (fname, "a"))) {
		fprintf (stderr, "ERROR: could not append field output.\n");
		return 0;
	}

	fwrite (field, sizeof(complex float), obs->count, fp);

	fclose (fp);

	return obs->count;
}

int getfields (char *fname, complex float *field, int len, float *nrm) {
	FILE *fp;
	complex float *fldptr;
	int nmeas, size[2], i, j;
	float lerr;

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: could not open field input.\n");
		return 0;
	}

	/* Read the number of recorded measurements. */
	fread (size, sizeof(int), 2, fp);

	nmeas = size[0] * size[1];

	/* The initial norm of the matrix is zero, if it is desired. */
	if (nrm) *nrm = 0;

	/* If the specified size doesn't match the recorded size, fail. */
	if (nmeas != len) {
		fprintf (stderr, "ERROR: recorded and specified counts do not agree.\n");
		return 0;
	}

	/* Read each column of the matrix and compute the norm. */
	for (fldptr = field, i = 0; i < size[1]; ++i) {
		fread (fldptr, sizeof(complex float), size[0], fp);
		if (nrm) {
			for (j = 0; j < size[0]; ++j) {
				lerr = cabs (fldptr[j]);
				*nrm += lerr * lerr;
			}
		}
		fldptr += size[0];
	}

	fclose (fp);

	return len;
}

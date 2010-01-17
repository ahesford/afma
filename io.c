#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "io.h"
#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"

#ifdef _FREEBSD
#define log2(a) (log(a) / log(2))
#endif

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

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
	int nmax;

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

	/* Compute the maximum FMM level for the desired finest level density. */
	fmaconf.maxlev = (int)ceil(log2(ceil((double)nmax / (double)fmaconf.bspbox)));

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
void getcontrast (char *fname, int *bslist, int nbs) {
	FILE *fp;
	int size[3], i, offset;
	long spos;

	/* Zero out the contrast in case nothing can be read. */
	memset (fmaconf.contrast, 0, nbs * sizeof(complex float));

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: unable to open %s.\n", fname);
		return;
	}

	/* Read the size of the contrast grid. */
	fread (size, sizeof(int), 3, fp);

	for (offset = 0, i = 0; i < nbs; ++i) {
		spos = (long)(bslist[i] - offset) * sizeof(complex float);
		offset = bslist[i] + 1;

		fseek (fp, spos, SEEK_CUR);
		fread (fmaconf.contrast + i, sizeof(complex float), 1, fp);
	}

	fclose (fp);
}

int prtcontrast (char *fname, complex float *currents) {
	int i, mpirank, size[3];
	FILE *fp;
	complex float *lct, *globct = NULL;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);

	lct = calloc (fmaconf.gnumbases, sizeof(complex float));

	if (!mpirank)
		globct = calloc (fmaconf.gnumbases, sizeof(complex float));

	for (i = 0; i < fmaconf.numbases; ++i)
		lct[fmaconf.bslist[i]] = currents[i];

	MPI_Reduce (lct, globct, 2 * fmaconf.gnumbases, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (mpirank) {
		free (lct);
		return fmaconf.gnumbases;
	}

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not open current output.\n");
		return 0;
	}

	size[0] = fmaconf.nx; size[1] = fmaconf.ny; size[2] = fmaconf.nz;
	fwrite (size, sizeof(int), 3, fp);
	fwrite (globct, sizeof(complex float), fmaconf.gnumbases, fp);

	fclose (fp);

	free (lct);
	free (globct);

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "io.h"
#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"
#include "integrate.h"

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

void gaqd_ (int *, double *, double *, double *, double *, int *);

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
	sscanf (buf, "%d", maxit);

	/* Read the regularization bounds and step. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %f", regparm, regparm + 1, regparm + 2);

	/* Read the DBIM tolerance. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", tol);

	fclose (fp);
}

/* Read the configuration file and set parameters. */
void getconfig (char *fname) {
	FILE *fp;
	int ierr;
	char buf[1024];

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
	sscanf (buf, "%f %f %f", fmaconf.min, fmaconf.min + 1, fmaconf.min + 2);

	/* Set the global number of bases, for easy reference. */
	fmaconf.gnumbases = fmaconf.nx * fmaconf.ny * fmaconf.nz;

	/* Read the upper corner of the domain. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %f", fmaconf.max, fmaconf.max + 1, fmaconf.max + 2);

	/* Compute the individual cell size. */
	fmaconf.cell[0] = (fmaconf.max[0] - fmaconf.min[0]) / fmaconf.nx;
	fmaconf.cell[1] = (fmaconf.max[1] - fmaconf.min[1]) / fmaconf.ny;
	fmaconf.cell[2] = (fmaconf.max[2] - fmaconf.min[2]) / fmaconf.nz;

	/* Set the wave number to 2 pi, since wavelength is the length unit. */
	fmaconf.k0 = 2 * M_PI;

	/* Read the MLFMA level configuration. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d %d", &(fmaconf.maxlev),
			&(fmaconf.toplev), &(fmaconf.sharedmax));

	/* Read the number of MLFMA buffer boxes. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.numbuffer));

	/* Read the MLFMA precision. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(fmaconf.precision));

	/* Read the MLFMA minimum box size. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(fmaconf.smallbox));

	/* Read the MLFMA interpolation order. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.interpord));

	/* Read the integration points for source and receiver. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d", &(fmaconf.nsrcint), &(fmaconf.nrcvint));

	/* Find the quadrature points for the source integration. */
	fmaconf.srcpts = malloc (2 * fmaconf.nsrcint * sizeof(double));
	fmaconf.srcwts = fmaconf.srcpts + fmaconf.nsrcint;
	gaqd_ (&(fmaconf.nsrcint), fmaconf.srcpts, fmaconf.srcwts, NULL, NULL, &ierr);
	for (ierr = 0; ierr < fmaconf.nsrcint; ++ierr)
		fmaconf.srcpts[ierr] = cos (fmaconf.srcpts[ierr]);

	/* Find the quadrature points for the receiver integration. */
	fmaconf.rcvpts = malloc (2 * fmaconf.nrcvint * sizeof(double));
	fmaconf.rcvwts = fmaconf.rcvpts + fmaconf.nrcvint;
	gaqd_ (&(fmaconf.nrcvint), fmaconf.rcvpts, fmaconf.rcvwts, NULL, NULL, &ierr);
	for (ierr = 0; ierr < fmaconf.nrcvint; ++ierr)
		fmaconf.rcvpts[ierr] = cos (fmaconf.rcvpts[ierr]);

	/* Read the iterative solver configuration. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d", &(solver.maxit), &(solver.restart));

	/* Read the iterative solver preconditioner setting. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(solver.precond));

	/* Read the iterative solver tolerance. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(solver.epscg));

	/* Read the source radius. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(srcmeas.radius));

	/* Read the source theta values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", srcmeas.trange, srcmeas.trange + 1, &(srcmeas.ntheta));

	/* Convert the degree values to radians. */
	srcmeas.trange[0] *= M_PI / 180;
	srcmeas.trange[1] *= M_PI / 180;

	/* Read the source phi values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", srcmeas.prange, srcmeas.prange + 1, &(srcmeas.nphi));

	/* Convert the degree values to radians. */
	srcmeas.prange[0] *= M_PI / 180;
	srcmeas.prange[1] *= M_PI / 180;

	/* Read the observer radius. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f", &(obsmeas.radius));

	/* Read the observer theta values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", obsmeas.trange, obsmeas.trange + 1, &(obsmeas.ntheta));

	/* Convert the degree values to radians. */
	obsmeas.trange[0] *= M_PI / 180;
	obsmeas.trange[1] *= M_PI / 180;

	/* Read the observer phi values. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%f %f %d", obsmeas.prange, obsmeas.prange + 1, &(obsmeas.nphi));

	/* Convert the degree values to radians. */
	obsmeas.prange[0] *= M_PI / 180;
	obsmeas.prange[1] *= M_PI / 180;

	fclose (fp);
}

/* Read a portion of the contrast file and store it. */
void getcontrast (char *fname, int *bslist, int nbs) {
	FILE *fp;
	int size[3], i, offset;
	long spos;

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

int prtfield (char *fname, measdesc *obs, complex float *field) {
	FILE *fp;
	int size[2];

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not open field output.\n");
		return 0;
	}

	size[0] = obs->nphi;
	size[1] = obs->ntheta;

	/* The locations are not stored. */
	fwrite (size, sizeof(int), 2, fp);
	fwrite (field, sizeof(complex float), obs->count, fp);

	fclose (fp);

	return obs->count;
}

int getfield (char *fname, complex float *field, int len) {
	FILE *fp;
	int nmeas, size[2];

	if (!(fp = fopen (fname, "r"))) {
		fprintf (stderr, "ERROR: could not open field input.\n");
		return 0;
	}

	/* Read the number of recorded measurements. */
	fread (size, sizeof(int), 2, fp);
	nmeas = size[0] * size[1];

	if (nmeas != len)
		fprintf (stderr, "ERROR: recorded and specified counts do not agree.\n");

	/* Read the observations. */
	fread (field, sizeof(complex float), MIN(len,nmeas), fp);
	fclose (fp);

	return len;
}

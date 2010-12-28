#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "precision.h"

#include "config.h"
#include "mlfma.h"
#include "itsolver.h"
#include "measure.h"

#include "util.h"

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
void getdbimcfg (char *fname, int *maxit, real *regparm, real *tol) {
	FILE *fp;
	char buf[1024];
	double rbuf[4];

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
	if (sscanf (buf, "%lf %lf %lf %lf", rbuf, rbuf + 1, rbuf + 2, rbuf + 3) < 4)
		regparm[3] = 1.0;
	else regparm[3] = (real) rbuf[3];
	/* Copy the remaining regularization parameters. */
	regparm[0] = (real) rbuf[0];
	regparm[1] = (real) rbuf[1];
	regparm[2] = (real) rbuf[2];

	/* Read the DBIM tolerance. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (sscanf (buf, "%lf %lf", rbuf, rbuf + 1) < 2) {
		tol[0] = (real) rbuf[0];
		tol[1] = tol[0] / 10.;
	} else {
		tol[0] = (real) rbuf[0];
		tol[1] = (real) rbuf[1];
	}

	fclose (fp);
}

/* Read the configuration file and set parameters. */
void getconfig (char *fname, solveparm *hislv, solveparm *loslv) {
	FILE *fp;
	char buf[1024];
	int nmax, nbox, i;
	double rbuf[4];

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
	sscanf (buf, "%lf %lf %lf %lf", rbuf, rbuf + 1, rbuf + 2, rbuf + 3);
	fmaconf.cen[0] = (real) rbuf[0];
	fmaconf.cen[1] = (real) rbuf[1];
	fmaconf.cen[2] = (real) rbuf[2];
	fmaconf.cell = (real) rbuf[3];

	fmaconf.min[0] = fmaconf.cen[0] - 0.5 * (real)(fmaconf.nx) * fmaconf.cell;
	fmaconf.min[1] = fmaconf.cen[1] - 0.5 * (real)(fmaconf.ny) * fmaconf.cell;
	fmaconf.min[2] = fmaconf.cen[2] - 0.5 * (real)(fmaconf.nz) * fmaconf.cell;

	fmaconf.cellvol = fmaconf.cell * fmaconf.cell * fmaconf.cell;

	/* Set the wave number to 2 pi, since wavelength is the length unit. */
	fmaconf.k0 = 2 * M_PI;

	/* Read the MLFMA level and fast translation configuration. If the fast
	 * O2I line isn't properly configured, just ignore fast O2I. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (sscanf (buf, "%d %d %d %d %d %d", &(fmaconf.bspbox),
				&(fmaconf.toplev), &(fmaconf.fo2itxlev),
				&(fmaconf.fo2ibclev), &(fmaconf.fo2iord),
				&(fmaconf.fo2iosr)) < 6)
		fmaconf.fo2itxlev = fmaconf.fo2ibclev = 
			fmaconf.fo2iord = fmaconf.fo2iosr = 0;

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

	/* Set defaults for fast, compressed O2I level designations. */
	/* Specifying a negative maximum expansion level means expand all. */
	if (fmaconf.fo2itxlev < 0) fmaconf.fo2itxlev = fmaconf.toplev;
	/* Specifying a negative minimum compression level means compress all. */
	if (fmaconf.fo2ibclev < 0) fmaconf.fo2ibclev = fmaconf.maxlev;

	/* Read the number of MLFMA buffer boxes. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.numbuffer));

	/* Read the MLFMA precision. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%lf", rbuf);
	fmaconf.precision = (real) rbuf[0];

	/* Read the MLFMA interpolation order. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d", &(fmaconf.interpord));

	/* Read the high-accuracy iterative solver configuration. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	sscanf (buf, "%d %d %lf", &(hislv->maxit), &(hislv->restart), rbuf);
	hislv->epscg = (real) rbuf[0];

	/* Read the low-accuracy iterative solver configuration, if desired. */
	skipcomments (fp);
	fgets (buf, 1024, fp);
	if (loslv) {
		sscanf (buf, "%d %d %lf", &(loslv->maxit), &(loslv->restart), rbuf);
		loslv->epscg = (real) rbuf[0];
	}

	fclose (fp);
}

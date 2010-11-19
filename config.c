#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

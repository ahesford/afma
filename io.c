#include <stdio.h>
#include <string.h>
#include <math.h>

#include "io.h"
#include "mlfma.h"
#include "itsolver.h"

void skipcomments (FILE *fp) {
	fpos_t loc;
	char buf[1024];

	do {
		fgetpos (fp, &loc);
		if (!fgets (buf, 1024, fp)) break;
	} while (strstr (buf, "#") != NULL);

	fsetpos (fp, &loc);
}

/* Read the configuration file and set parameters. */
void getconfig (char *fname) {
	FILE *fp;
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

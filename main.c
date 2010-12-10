#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <mpi.h>

#include "ScaleME.h"

#include "precision.h"

#include "itsolver.h"
#include "measure.h"
#include "direct.h"
#include "config.h"
#include "mlfma.h"
#include "util.h"
#include "io.h"

void usage (char *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-d] [-l #] [-r #] [-a #] [-n #] [-x] [-b] [-o <prefix>] -i <prefix>\n", name);
	fprintf (stderr, "\t-i: Specify input file prefix\n");
	fprintf (stderr, "\t-o: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "\t-d: Debug mode (prints induced field); specify twice to write after every restart\n");
	fprintf (stderr, "\t-b: Use BiCG-STAB instead of GMRES\n");
	fprintf (stderr, "\t-r: Specify the number of observation configurations\n");
	fprintf (stderr, "\t-a: Use ACA far-field transformations, or SVD when tolerance is negative\n");
	fprintf (stderr, "\t-l: Use loose GMRES with the specified number of augmented vectors\n");
	fprintf (stderr, "\t-n: Specify number of points for near-field integration\n");
	fprintf (stderr, "\t-x: Use singularity extraction for self terms\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist,
	     fname[1024], fldfmt[1024], guessfmt[1024];
	int mpirank, mpisize, i, j, k, nit, gsize[3], obscount = 1;
	int debug = 0, maxobs, useaca = 0, usebicg = 0, useloose = 0;
	int numsrcpts = 4, singex = 0;
	cplx *rhs, *sol, *field;
	double cputime, wtime;
	long nelt;
	real acatol = -1;
	augspace aug;

	measdesc *obsmeas, srcmeas;
	solveparm solver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	fprintf (stderr, "MPI Rank %d, pid %d\n", mpirank, getpid());

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:dbr:a:hl:xn:")) != -1) {
		switch (ch) {
		case 'i':
			inproj = optarg;
			break;
		case 'o':
			outproj = optarg;
			break;
		case 'd':
			++debug;
			break;
		case 'b':
			usebicg = 1;
			break;
		case 'r':
			obscount = strtol(optarg, NULL, 0);
			break;
		case 'a':
			acatol = strtod(optarg, NULL);
			useaca = 1;
			break;
		case 'l':
			useloose = strtol(optarg, NULL, 0);
			break;
		case 'x':
			singex = 1;
			break;
		case 'n':
			numsrcpts = strtol(optarg, NULL, 0);
			break;
		default:
			if (!mpirank) usage (arglist[0]);
			MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
		}
	}

	if (!inproj) {
		if (!mpirank) usage (arglist[0]);
		MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
	}

	if (!outproj) outproj = inproj;

	if (!mpirank) fprintf (stderr, "Reading configuration file.\n");

	/* Allocate the observation configurations. */
	obsmeas = calloc (obscount, sizeof(measdesc));

	/* Read the basic configuration. */
	sprintf (fname, "%s.input", inproj);
	getconfig (fname, &solver, NULL, &srcmeas, obsmeas, obscount);

	/* Convert the source range format to an explicit location list. */
	buildlocs (&srcmeas);
	/* Do the same for the observation locations. */
	for (i = 0; i < obscount; ++i) buildlocs (obsmeas + i);

	/* Initialize ScaleME and find the local basis set. */
	ScaleME_preconf (useaca);
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));

	nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	/* Allocate the RHS vector, which will also store the solution. */
	rhs = malloc (2 * nelt * sizeof(cplx));
	sol = rhs + nelt;
	/* Allocate the local portion of the contrast storage. */
	fmaconf.contrast = malloc (nelt * sizeof(cplx));
	/* Allocate the observation array. */
	for (i = 1, maxobs = obsmeas->count; i < obscount; ++i)
		maxobs = MAX(maxobs, obsmeas[i].count);
	field = malloc (maxobs * sizeof(cplx));

	/* Store the grid size for printing of field values. */
	gsize[0] = fmaconf.nx; gsize[1] = fmaconf.ny; gsize[2] = fmaconf.nz;

	wtime = MPI_Wtime();
	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the contrast for the local basis set. */
	sprintf (fname, "%s.contrast", inproj);
	/* Try both formats. */
	getctgrp (fmaconf.contrast, fname, gsize,
			fmaconf.bslist, fmaconf.numbases, fmaconf.bspbox);
	wtime = MPI_Wtime() - wtime;
	if (!mpirank) fprintf (stderr, "Wall time for contrast read: %0.6g\n", wtime);

	/* Precalculate some values for the FMM and direct interactions. */
	fmmprecalc (acatol, useaca);
	i = dirprecalc (numsrcpts, singex);
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	/* Build the root interpolation matrix for measurements. */
	for (i = 0; i < obscount; ++i)
		ScaleME_buildRootInterpMat (obsmeas[i].imat, fmaconf.interpord,
				obsmeas[i].ntheta, obsmeas[i].nphi,
				obsmeas[i].trange, obsmeas[i].prange);

	/* Find the width of the integer label in the field name. */
	i = (int)ceil(log10(srcmeas.count));
	if (obscount < 2)
		sprintf (fldfmt, "%%s.tx%%0%dd.field", i);
	else {
		j = (int)ceil(log10(obscount));
		sprintf (fldfmt, "%%s.tx%%0%dd.rx%%0%dd.field", i, j);
	}

	sprintf (guessfmt, "%%s.tx%%0%dd.%%s", i);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	/* Ensure each process is waiting at the start of the loop. */
	MPI_Barrier (MPI_COMM_WORLD);

	/* Initialize the loose GMRES buffer. */
	if (useloose > 0) {
		aug.start = -1;
		aug.nmax = useloose;
		aug.ntot = 0;
		aug.z = malloc(2 * aug.nmax * nelt * sizeof(cplx));
		aug.az = aug.z + aug.nmax * nelt;
	}

	for (i = 0; i < srcmeas.count; ++i) {

		if (!mpirank)
			fprintf (stderr, "Running simulation for source %d.\n", i + 1);

		/* Attempt to read the pre-computed RHS from a file.
		 * If this fails, build the RHS directly. */
		sprintf (fname, guessfmt, inproj, i, "rhs");
		if (!getctgrp (rhs, fname, gsize, fmaconf.bslist,
					fmaconf.numbases, fmaconf.bspbox)) {
			if (!mpirank) fprintf (stderr, "Building RHS.\n");
			buildrhs (rhs, srcmeas.locations + 3 * i);
		}

		/* Attempt to read an initial first guess from a file. */
		sprintf (fname, guessfmt, inproj, i, "guess");
		k = getctgrp (sol, fname, gsize, fmaconf.bslist,
				fmaconf.numbases, fmaconf.bspbox);

		/* Restart if the true residual is not sufficiently low. */
		for (j = 0, nit = 1; j < solver.restart && nit > 0; ++j) {
			cputime = (double)clock() / CLOCKS_PER_SEC;
			wtime = MPI_Wtime();

			if (usebicg) nit = bicgstab (rhs, sol, k || j,
					solver.maxit, solver.epscg, 0);
			else nit = gmres (rhs, sol, k || j, solver.maxit,
					solver.epscg, 0, useloose > 0 ? &aug : NULL);

			cputime = (double)clock() / CLOCKS_PER_SEC - cputime;
			wtime = MPI_Wtime() - wtime;

			if (!mpirank) {
				fprintf (stderr, "CPU time for solution: %0.6g\n", cputime);
				fprintf (stderr, "Wall time for solution: %0.6g\n", wtime);
			}

			/* Update the internal field before every restart, but
			 * only when the debug flag has been specified twice. */
			if (debug > 1) {
				wtime = MPI_Wtime();
				sprintf (fname, guessfmt, outproj, i, "solution");
				prtctgrp (fname, sol, gsize, fmaconf.bslist,
						fmaconf.numbases, fmaconf.bspbox);
				wtime = MPI_Wtime() - wtime;
				if (!mpirank) fprintf (stderr, "Wall time for field write: %0.6g\n", wtime);
			}
		}

		/* Update the field after all iterations, but only with a
		 * single debug flag. Otherwise, it has already been written. */
		if (debug == 1) {
			wtime = MPI_Wtime();
			sprintf (fname, guessfmt, outproj, i, "solution");
			prtctgrp (fname, sol, gsize, fmaconf.bslist,
					fmaconf.numbases, fmaconf.bspbox);
			wtime = MPI_Wtime() - wtime;
			if (!mpirank) fprintf (stderr, "Wall time for field write: %0.6g\n", wtime);
		}

		/* Convert total field into contrast current. */
		for (j = 0; j < nelt; ++j) sol[j] *= fmaconf.contrast[j];

		/* Compute and write out all observation fields. */
		for (k = 0; k < obscount; ++k) {
			farfield (sol, obsmeas + k, field);

			/* Append the field for the current transmitter. */
			if (!mpirank) {
				sprintf (fname, fldfmt,  outproj, i, k);
				writefld (fname, obsmeas[k].nphi, obsmeas[k].ntheta, field);
			}
		}
	}

	ScaleME_finalizeParHostFMA ();

	freedircache ();
	delmeas (&srcmeas);
	for (i = 0; i < obscount; ++i) delmeas (obsmeas + i);
	free (obsmeas);

	if (useloose > 0) free (aug.z);

	free (rhs);
	free (field);
	free (fmaconf.contrast);
	free (fmaconf.radpats);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

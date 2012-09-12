#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_max_threads () { return 1; }
#endif /* _OPENMP */

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
	fprintf (stderr, "Usage: %s [-d] [-l #] [-a #] [-n #] [-b] [-f x,y,z,a]\n"
			 "       [-o <prefix>] -s <src> -r <obs> -i <prefix>\n", name);
	fprintf (stderr, "  -i: Specify input file prefix\n");
	fprintf (stderr, "  -o: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "  -d: Debug mode (prints induced field); specify twice to write after every restart\n");
	fprintf (stderr, "  -b: Use BiCG-STAB instead of GMRES\n");
	fprintf (stderr, "  -a: Use ACA far-field transformations, or SVD when tolerance is negative\n");
	fprintf (stderr, "  -l: Use loose GMRES with the specified number of augmented vectors\n");
	fprintf (stderr, "  -n: Specify number of points for near-field integration\n");
	fprintf (stderr, "  -s: Specify the source location or range\n");
	fprintf (stderr, "  -r: Specify the observation range\n");
	fprintf (stderr, "  -f: Specify a focal axis x,y,z and width a for the incident field\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024],
	     fldfmt[1024], guessfmt[1024], *srcspec = NULL, *obspec = NULL;
	int mpirank, mpisize, i, j, k, nit, gsize[3];
	int debug = 0, maxobs, useaca = 0, usebicg = 0, useloose = 0, usedir = 0;
	int numsrcpts = 5;
	cplx *rhs, *sol, *inc, *field;
	double cputime, wtime;
	long nelt;
	real acatol = -1;
	real dir[4];
	augspace aug;

	measdesc obsmeas, srcmeas;
	solveparm solver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	fprintf (stderr, "MPI Rank %d, pid %d, %d threads\n",
			mpirank, getpid(), omp_get_max_threads());

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:dba:hl:n:s:r:f:")) != -1) {
		switch (ch) {
		case 'i':
			inproj = optarg;
			break;
		case 'o':
			outproj = optarg;
			break;
		case 'd':
			debug = 1;
			break;
		case 'b':
			usebicg = 1;
			break;
		case 'a':
			acatol = strtod(optarg, NULL);
			useaca = 1;
			break;
		case 'l':
			useloose = strtol(optarg, NULL, 0);
			break;
		case 'n':
			numsrcpts = strtol(optarg, NULL, 0);
			break;
		case 'r':
			obspec = optarg;
			break;
		case 's':
			srcspec = optarg;
			break;
		case 'f':
			dir[0] = strtod(strtok(optarg, ","), NULL);
			dir[1] = strtod(strtok(NULL, ","), NULL);
			dir[2] = strtod(strtok(NULL, ","), NULL);
			dir[3] = strtod(strtok(NULL, ","), NULL);
			usedir = 1;
			break;
		default:
			if (!mpirank) usage (arglist[0]);
			MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
		}
	}

	/* The project name must be specified. */
	if (!inproj || !srcspec || !obspec) {
		if (!mpirank) usage (arglist[0]);
		MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
	}

	if (!outproj) outproj = inproj;

	if (!mpirank) fprintf (stderr, "Reading configuration file.\n");

	/* Read the basic configuration. */
	sprintf (fname, "%s.input", inproj);
	getconfig (fname, &solver, NULL);

	/* Build the source and observer location specifiers. */
	buildsrc (&srcmeas, srcspec);
	buildobs (&obsmeas, obspec);

	/* Initialize ScaleME and find the local basis set. */
	ScaleME_preconf (useaca);
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));

	nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	/* Allocate the RHS vector, solution and incident field. */
	rhs = malloc (3 * nelt * sizeof(cplx));
	sol = rhs + nelt;
	inc = sol + nelt;
	/* Allocate the local portion of the contrast storage. */
	fmaconf.contrast = malloc (nelt * sizeof(cplx));
	/* Allocate the observation array. */
	field = malloc (obsmeas.count * sizeof(cplx));

	/* Store the grid size for printing of field values. */
	gsize[0] = fmaconf.nx; gsize[1] = fmaconf.ny; gsize[2] = fmaconf.nz;

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the contrast for the local basis set. */
	sprintf (fname, "%s.contrast", inproj);
	/* Try both formats. */
	getctgrp (fmaconf.contrast, fname, gsize,
			fmaconf.bslist, fmaconf.numbases, fmaconf.bspbox);

	/* First build the integration rules that will be used. */
	bldintrules (numsrcpts, 0);
	/* Precalculate some values for the FMM and direct interactions. */
	fmmprecalc (acatol, useaca);
	i = dirprecalc (numsrcpts);
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	/* Build the root interpolation matrix for measurements. */
	ScaleME_buildRootInterpMat (obsmeas.imat, fmaconf.interpord,
			obsmeas.ntheta, obsmeas.nphi, obsmeas.trange, obsmeas.prange);

	/* Find the width of the integer label in the field name. */
	i = (int)ceil(log10(srcmeas.count));
	sprintf (fldfmt, "%%s.tx%%0%dd.farfld", i);
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
			/* Integrate the incident field over each cell. */
			buildrhs (inc, srcmeas.locations + 3 * i, srcmeas.plane, usedir ? dir : NULL);
			/* Convert the integrated field to the average. */
#pragma omp parallel for default(shared) private(j)
			for (j = 0; j < nelt; ++j) inc[j] /= fmaconf.cellvol;
			/* Apply the scattering integral to the incident field. */
			matvec(rhs, inc, sol, 0);
		} else {
			/* If the RHS was read, it is the incident field.
			 * Set the incident field in inc to zero. */
#pragma omp parallel for default(shared) private(j)
			for (j = 0; j < nelt; ++j) inc[j] = 0.0;
		}

		/* Attempt to read an initial first guess from a file. */
		sprintf (fname, guessfmt, inproj, i, "guess");
		k = getctgrp (sol, fname, gsize, fmaconf.bslist,
				fmaconf.numbases, fmaconf.bspbox);

		/* Restart if the true residual is not sufficiently low. */
		for (j = 0, nit = 1; j < solver.restart && nit > 0; ++j) {
			/* Reset the debug tripwire to ensure the final output
			 * is written if intermediate writes were performed. */
			debug = debug ? 1 : 0;

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

			/* Update the total field before every restart, but
			 * only if the tripwire file exists. */
			sprintf (fname, guessfmt, inproj, i, "volwrite");
			if (!access (fname, F_OK)) {
				/* The master process should unlink the tripwire. */
				MPI_Barrier (MPI_COMM_WORLD);
				if (!mpirank) unlink (fname);
				sprintf (fname, guessfmt, outproj, i, "scatfld");
				prtctgrp (fname, sol, gsize, fmaconf.bslist,
						fmaconf.numbases, fmaconf.bspbox);
				/* If restarts have finished, avoid an
				 * immediate rewrite of this field. */
				debug = debug ? 2 : 0;
			}
		}

		/* Write the scattered field if desired and not just written. */
		if (debug == 1) {
			sprintf (fname, guessfmt, outproj, i, "scatfld");
			prtctgrp (fname, sol, gsize, fmaconf.bslist,
					fmaconf.numbases, fmaconf.bspbox);
		}

		/* Convert the scattered field to the total field. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) sol[j] += inc[j];

		/* Write the total field if desired and not just written. */
		if (debug > 0) {
			sprintf (fname, guessfmt, outproj, i, "totfld");
			prtctgrp (fname, sol, gsize, fmaconf.bslist,
					fmaconf.numbases, fmaconf.bspbox);
		}

		/* Convert total field into contrast current. */
#pragma omp parallel for default(shared) private(j)
		for (j = 0; j < nelt; ++j) sol[j] *= fmaconf.contrast[j];

		/* Compute and write out all observation fields. */
		farfield (sol, &obsmeas, field);

		/* Append the field for the current transmitter. */
		if (!mpirank) {
			sprintf (fname, fldfmt,  outproj, i);
			writefld (fname, obsmeas.nphi, obsmeas.ntheta, field);
		}
	}

	ScaleME_finalizeParHostFMA ();

	freedircache ();
	delmeas (&srcmeas);
	delmeas (&obsmeas);
	delintrules ();

	if (useloose > 0) free (aug.z);

	free (rhs);
	free (field);
	free (fmaconf.contrast);
	free (fmaconf.radpats);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <mpi.h>

#include "ScaleME.h"

#include "itsolver.h"
#include "measure.h"
#include "excite.h"
#include "direct.h"
#include "mlfma.h"
#include "io.h"

void usage (char *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-d] [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "\t-d: Debug mode (prints RHS and induced currents)\n");
	fprintf (stderr, "\t-g: Use GMRES instead of BiCG-STAB\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, j, nit, gmr = 0, gsize[3];
	complex float *rhs, *sol, *field;
	double cputime, wtime;
	int debug = 0;
	long nelt;

	measdesc obsmeas, srcmeas;
	solveparm solver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	fprintf (stderr, "MPI Rank %d, pid %d\n", mpirank, getpid());

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:dg")) != -1) {
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
		case 'g':
			gmr = 1;
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

	/* Read the basic configuration. */
	sprintf (fname, "%s.input", inproj);
	getconfig (fname, &solver, NULL, &srcmeas, &obsmeas);

	/* Convert the source range format to an explicit location list. */
	buildlocs (&srcmeas);
	/* Do the same for the observation locations. */
	buildlocs (&obsmeas);

	/* Initialize ScaleME and find the local basis set. */
	ScaleME_preconf ();
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));

	nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;
	/* Allocate the RHS vector, which will also store the solution. */
	rhs = malloc (2 * nelt * sizeof(complex float));
	sol = rhs + nelt;
	/* Allocate the local portion of the contrast storage. */
	fmaconf.contrast = malloc (nelt * sizeof(complex float));
	/* Allocate the observation array. */
	field = malloc (obsmeas.count * sizeof(complex float));

	/* Store the grid size for printing of field values. */
	gsize[0] = fmaconf.nx; gsize[1] = fmaconf.ny; gsize[2] = fmaconf.nz;

	wtime = MPI_Wtime();
	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the contrast for the local basis set. */
	sprintf (fname, "%s.contrast", inproj);
	getcontrast (fmaconf.contrast, fname, gsize,
			fmaconf.bslist, fmaconf.numbases, fmaconf.bspbox);
	wtime = MPI_Wtime() - wtime;
	if (!mpirank) fprintf (stderr, "Wall time for contrast read: %0.6g\n", wtime);

	/* Precalculate some values for the FMM and direct interactions. */
	fmmprecalc ();
	i = dirprecalc ();
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");
	
	for (i = 0; i < srcmeas.count; ++i) {
		if (!mpirank)
			fprintf (stderr, "Running simulation for source %d.\n", i + 1);
		/* Build the right-hand side for the specified location. Use
		 * point sources, rather than plane waves, for excitation. */
		buildrhs (rhs, srcmeas.locations + 3 * i);

		/* Initial first guess is zero. */
		memset (sol, 0, nelt * sizeof(complex float));

		cputime = (double)clock() / CLOCKS_PER_SEC;
		wtime = MPI_Wtime();
		/* Use GMRES if requested, otherwise use BiCG-STAB. */
		if (gmr) cgmres (rhs, sol, 0, &solver);
		else {
			/* BiCG-STAB restarts if the true residual differs from
			 * the predicted residual. */
			for (j = 0, nit = 1; j < solver.restart && nit > 0; ++j)
				nit = bicgstab (rhs, sol, j, solver.maxit, solver.epscg);
		}
		cputime = (double)clock() / CLOCKS_PER_SEC - cputime;
		wtime = MPI_Wtime() - wtime;

		if (!mpirank) {
			fprintf (stderr, "CPU time for solution: %0.6g\n", cputime);
			fprintf (stderr, "Wall time for solution: %0.6g\n", wtime);
		}

		if (debug) {
			wtime = MPI_Wtime();
			sprintf (fname, "%s.%d.currents", outproj, i);
			prtcontrast (fname, sol, gsize, fmaconf.bslist,
					fmaconf.numbases, fmaconf.bspbox);
			wtime = MPI_Wtime() - wtime;
			if (!mpirank) fprintf (stderr, "Wall time for field write: %0.6g\n", wtime);
		}

		/* Convert total field into contrast current. */
		for (j = 0; j < nelt; ++j)
			sol[j] *= fmaconf.contrast[j];

		farfield (sol, &obsmeas, field);

		/* Append the field for the current transmitter. */
		if (!mpirank) {
			sprintf (fname, "%s.%d.field", outproj, i);
			writefld (fname, &obsmeas, field);
		}
	}

	ScaleME_finalizeParHostFMA ();

	freedircache ();

	free (rhs);
	free (field);
	free (fmaconf.contrast);
	free (fmaconf.radpats);
	free (srcmeas.locations);
	free (obsmeas.locations);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

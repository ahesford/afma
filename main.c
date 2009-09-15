#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <mpi.h>

#include "ScaleME.h"

#include "itsolver.h"
#include "measure.h"
#include "parfmm.h"
#include "excite.h"
#include "mlfma.h"
#include "io.h"

void usage (char *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-d] [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "\t-d: Debug mode (prints RHS and induced currents)\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, j;
	complex float *rhs, *field;
	clock_t tstart, tend;
	double cputime;
	int debug = 0;

	measdesc obsmeas, srcmeas;
	solveparm solver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	fprintf (stderr, "MPI Rank %d, pid %d\n", mpirank, getpid());

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:d")) != -1) {
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
	/* Allocate the RHS vector, which will also store the solution. */
	rhs = malloc (fmaconf.numbases * sizeof(complex float));
	/* Allocate the local portion of the contrast storage. */
	fmaconf.contrast = malloc (fmaconf.numbases * sizeof(complex float));
	/* Allocate the observation array. */
	field = malloc (obsmeas.count * sizeof(complex float));

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the contrast for the local basis set. */
	sprintf (fname, "%s.contrast", inproj);
	getcontrast (fname, fmaconf.bslist, fmaconf.numbases);

	i = preimpedance ();
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	sprintf (fname, "%s.field", outproj);
	prtfldhdr (fname, &srcmeas, &obsmeas);

	for (i = 0; i < srcmeas.count; ++i) {
		if (!mpirank)
			fprintf (stderr, "Running simulation for source %d.\n", i + 1);
		/* Build the right-hand side for the specified location. Use
		 * point sources, rather than plane waves, for excitation. */
		buildrhs (rhs, srcmeas.locations + 3 * i);

		if (debug) {
			sprintf (fname, "%s.%d.rhs", outproj, i);
			prtcontrast (fname, rhs);
		}

		tstart = clock ();
		/* Run the iterative solver. The solution is stored in the RHS. */
		cgmres (rhs, rhs, 0, &solver);
		tend = clock ();
		cputime = (double) (tend - tstart) / CLOCKS_PER_SEC;

		if (!mpirank) fprintf (stderr, "CPU time for CGMRES: %f\n", cputime);

		if (debug) {
			sprintf (fname, "%s.%d.currents", outproj, i);
			prtcontrast (fname, rhs);
		}

		/* Convert total field into contrast current. */
		for (j = 0; j < fmaconf.numbases; ++j)
			rhs[j] *= fmaconf.contrast[j];

		farfield (rhs, &obsmeas, field);

		/* Append the field for the current transmitter. */
		if (!mpirank) {
			sprintf (fname, "%s.field", outproj);
			appendfld (fname, &obsmeas, field);
		}
	}

	ScaleME_finalizeParHostFMA ();

	free (rhs);
	free (field);
	free (fmaconf.contrast);
	free (srcmeas.locations);
	free (obsmeas.locations);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

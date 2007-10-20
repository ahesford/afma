#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>

#include <mpi.h>

#include <ScaleME.h> // ScaleME-provided header
#include <Complex.h> // ScaleME complex header

#include "itsolver.h"
#include "measure.h"
#include "scaleme.h"
#include "excite.h"
#include "mlfma.h"
#include "io.h"

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i;
	complex float *rhs;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:")) != -1) {
		switch (ch) {
		case 'i':
			inproj = optarg;
			break;
		case 'o':
			outproj = optarg;
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
	getconfig (fname);

	/* Convert the source range format to an explicit location list. */
	buildlocs (&srcmeas);

	/* Initialize ScaleME and find the local basis set. */
	ScaleME_preconf ();
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));
	/* Allocate the RHS vector, which will also store the solution. */
	rhs = malloc (fmaconf.numbases * sizeof(complex float));

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the contrast for the local basis set. */
	sprintf (fname, "%s.contrast", inproj);
	getcontrast (fname, fmaconf.bslist, fmaconf.numbases);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	for (i = 0; i < srcmeas.count; ++i) {
		if (!mpirank)
			fprintf (stderr, "Running simulation for source %d.\n", i + 1);
		/* Build the right-hand side for the specified location. Use
		 * point sources, rather than plane waves, for excitation. */
		buildrhs (rhs, srcmeas.locations + 3 * i, 0);

		MPI_Barrier (MPI_COMM_WORLD);

		/* Run the iterative solver. The solution is stored in the RHS. */
		cgmres (rhs, rhs);

		MPI_Barrier (MPI_COMM_WORLD);
	}

	ScaleME_finalizeParHostFMA ();

	free (rhs);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

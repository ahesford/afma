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

void usage (char *);
int printcrt (char *, complex float *);
int printfield (char *, measdesc *, complex float *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
}

int printcrt (char *fname, complex float *currents) {
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

	if (mpirank) return fmaconf.gnumbases;

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not open current output.\n");
		return 0;
	}

	size[0] = fmaconf.nx; size[1] = fmaconf.ny; size[2] = fmaconf.nz;
	fwrite (size, sizeof(int), 3, fp);
	fwrite (globct, sizeof(complex float), fmaconf.gnumbases, fp);

	fclose (fp);

	return fmaconf.gnumbases;
}

int printfield (char *fname, measdesc *obs, complex float *field) {
	FILE *fp;

	if (!(fp = fopen (fname, "w"))) {
		fprintf (stderr, "ERROR: could not open current output.\n");
		return 0;
	}

	fwrite (&(obs->count), sizeof(int), 1, fp);
	fwrite (obs->locations, sizeof(float), 3 * obs->count, fp);
	fwrite (field, sizeof(complex float), obs->count, fp);

	fclose (fp);

	return obs->count;
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, j;
	complex float *rhs, *field;

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

		sprintf (fname, "%s.%d.rhs", outproj, i);
		printcrt (fname, rhs);

		MPI_Barrier (MPI_COMM_WORLD);

		/* Run the iterative solver. The solution is stored in the RHS. */
		cgmres (rhs, rhs);

		/* Convert total field into contrast current. */
		for (j = 0; j < fmaconf.numbases; ++j)
			rhs[i] *= fmaconf.contrast[i];

		sprintf (fname, "%s.%d.currents", outproj, i);
		printcrt (fname, rhs);

		MPI_Barrier (MPI_COMM_WORLD);

		if (obsmeas.radius < 10) directfield (rhs, &obsmeas, field);
		else farfield (rhs, &obsmeas, field);

		if (!mpirank) {
			sprintf (fname, "%s.%d.field", outproj, i);
			printfield (fname, &obsmeas, field);
		}

		MPI_Barrier (MPI_COMM_WORLD);
	}

	ScaleME_finalizeParHostFMA ();

	free (rhs);
	free (field);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

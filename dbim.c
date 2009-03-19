#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include <ScaleME.h> // ScaleME-provided header
#include <Complex.h> // ScaleME complex header

#include "itsolver.h"
#include "measure.h"
#include "frechet.h"
#include "scaleme.h"
#include "excite.h"
#include "mlfma.h"
#include "io.h"
#include "cg.h"

void usage (char *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, j, k, nmeas, dbimit;
	complex float *rhs, *rn, *currents, *field, *error, *errptr, *fldptr;
	float errnorm, tolerance, lerr, regparm[3], cgnorm, erninc;

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

	/* Read the DBIM-specific configuration. */
	sprintf (fname, "%s.dbimin", inproj);
	getdbimcfg (fname, &dbimit, regparm, &tolerance);

	/* Convert the source range format to an explicit location list. */
	buildlocs (&srcmeas);
	/* Do the same for the observation locations. */
	buildlocs (&obsmeas);
	/* The total number of measurements. */
	nmeas = srcmeas.count * obsmeas.count;

	/* Initialize ScaleME and find the local basis set. */
	ScaleME_preconf ();
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));

	/* Allocate the RHS vector, residual vector and contrast. */
	rhs = malloc (4 * fmaconf.numbases * sizeof(complex float));
	fmaconf.contrast = rhs + fmaconf.numbases;
	rn = fmaconf.contrast + fmaconf.numbases;
	currents = rn + fmaconf.numbases;

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the guess contrast for the local basis set. */
	sprintf (fname, "%s.guess", inproj);
	getcontrast (fname, fmaconf.bslist, fmaconf.numbases);

	i = preimpedance ();
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	/* Allocate the observation array. */
	field = malloc (2 * nmeas * sizeof(complex float));
	error = field + nmeas;

	/* Read the measurements and compute their norm. */
	sprintf (fname, "%s.field", inproj);
	getfields (fname, field, nmeas, &erninc);

	bldfrechbuf (fmaconf.numbases);

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	for (i = 0; i < dbimit; ++i) {
		if (!mpirank) fprintf (stderr, "DBIM iteration %d.\n", i);

		errptr = error;
		fldptr = field;
		errnorm = 0;
		memset (rn, 0, fmaconf.numbases * sizeof(complex float));
		for (j = 0; j < srcmeas.count; ++j) {
			/* Build the right-hand side for the specified location. Use
			 * point sources, rather than plane waves, for excitation. */
			buildrhs (rhs, srcmeas.locations + 3 * j);

			MPI_Barrier (MPI_COMM_WORLD);
			/* Run the iterative solver. The solution is stored in the RHS. */
			cgmres (rhs, rhs, 1); 

			/* Convert total field into contrast current. */
			for (k = 0; k < fmaconf.numbases; ++k)
				currents[k] = rhs[k] * fmaconf.contrast[k]; 
			
			MPI_Barrier (MPI_COMM_WORLD);

			/* Evaluate the scattered field. */
			farfield (currents, &obsmeas, errptr);
			MPI_Bcast (errptr, 2 * obsmeas.count, MPI_FLOAT, 0, MPI_COMM_WORLD);

			/* Compute the error vector. */
			for (k = 0; k < obsmeas.count; ++k) {
				errptr[k] = fldptr[k] - errptr[k];
				lerr = cabs(errptr[k]);
				errnorm += lerr * lerr;
			}

			/* Evaluate the adjoint Frechet derivative. */
			frechadj (errptr, rhs, rn);

			errptr += obsmeas.count;
			fldptr += obsmeas.count;
		}

		errnorm = sqrt (errnorm  / erninc);

		if (!mpirank)
			fprintf (stderr, "DBIM relative error: %f, iteration %d.\n", errnorm, i);
		if (errnorm < tolerance) break;

		/* Solve the system with CG for least-squares. */
		cgnorm = cgls (rn, currents, regparm[0]);

		if (!mpirank)
			fprintf (stderr, "CG relative error: %f, iteration %d.\n", cgnorm, i);

		/* Update the background. */
		for (j = 0; j < fmaconf.numbases; ++j)
			fmaconf.contrast[j] += currents[j]; 
		
		sprintf (fname, "%s.inverse.%03d", outproj, i);
		prtcontrast (fname, fmaconf.contrast);

		/* Scale the DBIM regularization parameter, if appropriate. */
		if (regparm[0] > regparm[1]) regparm[0] *= regparm[2];
	}

	ScaleME_finalizeParHostFMA ();

	free (rhs);
	free (field);
	delfrechbuf ();

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

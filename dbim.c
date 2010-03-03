#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "ScaleME.h"

#include "itsolver.h"
#include "measure.h"
#include "frechet.h"
#include "excite.h"
#include "direct.h"
#include "mlfma.h"
#include "io.h"
#include "cg.h"

extern complex float *rcvgrf;

void usage (char *);
float dbimerr (complex float *, complex float *, complex float *,
		solveparm *, solveparm *, measdesc *, measdesc *);

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-l] [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "\t-l: Disable default leapfrog (Kaczmarz-like) iterations\n");
}

float dbimerr (complex float *error, complex float *rn, complex float *field,
	solveparm *hislv, solveparm *loslv, measdesc *src, measdesc *obs) {
	complex float *rhs, *crt, *err, *fldptr;
	float errnorm = 0, lerr, errd = 0;
	int j, k;

	if (!error) err = malloc (obs->count * sizeof(complex float));
	else err = error;

	rhs = malloc (2 * fmaconf.numbases * sizeof(complex float));
	crt = rhs + fmaconf.numbases;

	if (rn) memset (rn, 0, fmaconf.numbases * sizeof(complex float));
	
	for (j = 0, fldptr = field; j < src->count; ++j, fldptr += obs->count) {
		/* Build the right-hand side for the specified location. Use
		 * point sources, rather than plane waves, for excitation. */
		buildrhs (rhs, src->locations + 3 * j);
		
		MPI_Barrier (MPI_COMM_WORLD);
		/* Run the iterative solver. The solution is stored in the RHS. */
		cgmres (rhs, rhs, 1, hislv);
		
		/* Convert total field into contrast current. */
		for (k = 0; k < fmaconf.numbases; ++k)
			crt[k] = rhs[k] * fmaconf.contrast[k];
		
		MPI_Barrier (MPI_COMM_WORLD);
		
		/* Evaluate the scattered field. */
		farfield (crt, obs, err);
		
		/* Compute the error vector. */
		for (k = 0; k < obs->count; ++k) {
			err[k] = fldptr[k] - err[k];
			lerr = cabs(err[k]);
			errnorm += lerr * lerr;
			lerr = cabs(fldptr[k]);
			errd += lerr * lerr;
		}
		
		/* Evaluate the adjoint Frechet derivative. */
		if (rn) frechadj (err, rhs, rn, obs, loslv);
		/* Save this error column, if storage is available. */
		if (error) err += obs->count;
	}

	free (rhs);
	if (!error && err) free (err);

	return sqrt(errnorm / errd);
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, j, nmeas, dbimit[2], q, lfrog = 1, gsize[3];
	complex float *rn, *crt, *field, *fldptr, *error;
	float errnorm = 0, tolerance[2], regparm[4], cgnorm, erninc;
	solveparm hislv, loslv;
	measdesc obsmeas, srcmeas, ssrc;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:l")) != -1) {
		switch (ch) {
		case 'i':
			inproj = optarg;
			break;
		case 'o':
			outproj = optarg;
			break;
		case 'l':
			lfrog = 0;
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
	getconfig (fname, &hislv, &loslv, &srcmeas, &obsmeas);

	/* Read the DBIM-specific configuration. */
	sprintf (fname, "%s.dbimin", inproj);
	getdbimcfg (fname, dbimit, regparm, tolerance);

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
	fmaconf.contrast = malloc (3 * fmaconf.numbases * sizeof(complex float));
	rn = fmaconf.contrast + fmaconf.numbases;
	crt = rn + fmaconf.numbases;

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the guess contrast for the local basis set. */
	sprintf (fname, "%s.guess", inproj);
	getcontrast (fmaconf.contrast, fname, fmaconf.bslist, fmaconf.numbases);

	/* Store the 3-D grid size for writing the contrast values. */
	gsize[0] = fmaconf.nx; gsize[1] = fmaconf.ny; gsize[2] = fmaconf.nz;

	/* Precalculate some values for the FMM and direct interactions. */
	fmmprecalc ();
	i = dirprecalc ();
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	/* Allocate the observation array. */
	field = malloc (nmeas * sizeof(complex float));

	/* Read the measurements and compute their norm. */
	sprintf (fname, "%s.field", inproj);
	getfields (fname, field, nmeas, &erninc);

	bldfrechbuf (fmaconf.numbases, &obsmeas);

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	/* One source per iteration. */
	ssrc.count = 1;

	/* The error storage vector for a single transmitter. */
	error = malloc (obsmeas.count * sizeof(complex float));

	/* Start a two-pass DBIM with low tolerances first. */
	if (!mpirank) fprintf (stderr, "First DBIM pass (%d iterations)\n", dbimit[0]);
	for (i = 0; i < dbimit[0]; ++i) {
		if (lfrog) {
			/* Use the leapfrogging (Kaczmarz-like) method. */
			for (q = 0, fldptr = field; q < srcmeas.count; ++q, fldptr += obsmeas.count) {
				ssrc.locations = srcmeas.locations + 3 * q;
				errnorm = dbimerr (error, NULL, fldptr, &hislv, &loslv, &ssrc, &obsmeas);
				
				/* Solve the system with CG for minimum norm. */
				cgnorm = cgmn (error, crt, &loslv, &ssrc, &obsmeas, regparm[0]);
				
				if (!mpirank)
					fprintf (stderr, "DBIM: %g, CG: %g (%d/%d).\n", errnorm, cgnorm, i, q);
				
				/* Update the background. */
				for (j = 0; j < fmaconf.numbases; ++j)
					fmaconf.contrast[j] += crt[j];
				
				sprintf (fname, "%s.inverse.t%03d", outproj, q);
				prtcontrast (fname, fmaconf.contrast, 
						fmaconf.bslist, fmaconf.numbases, gsize);
			}
			
			if (!mpirank) fprintf (stderr, "Reassess DBIM error.\n");
			errnorm = dbimerr (NULL, NULL, field, &hislv, &loslv, &srcmeas, &obsmeas);
		} else {
			errnorm = dbimerr (NULL, rn, field, &hislv, &loslv, &srcmeas, &obsmeas);
			/* Solve the system with CG for least squares. */
			cgnorm = cgls (rn, crt, &loslv, &srcmeas, &obsmeas, regparm[0]);
			
			/* Update the background. */
			for (j = 0; j < fmaconf.numbases; ++j)
				fmaconf.contrast[j] += crt[j];
		}
		
		sprintf (fname, "%s.inverse.%03d", outproj, i);
		prtcontrast (fname, fmaconf.contrast,
				fmaconf.bslist, fmaconf.numbases, gsize);
		
		if (!mpirank)
			fprintf (stderr, "DBIM relative error: %g, iteration %d.\n", errnorm, i);
		if (errnorm < tolerance[0]) break;
		
		/* Scale the DBIM regularization parameter, if appropriate. */
		if (!((i + 1) % (int)(regparm[3])) && regparm[0] > regparm[1]) {
			regparm[0] *= regparm[2];
			if (!mpirank)
				fprintf (stderr, "Regularization parameter: %g\n", regparm[0]);
		}
	}

	free (error);
	error = NULL;

	/* More iterations at high accuracy, without leapfrogging. */
	if (i < dbimit[0]) ++i;

	/* If the number of measurements is less than the number of pixels,
	 * use the minimum-norm conjugate gradient. */
	if (nmeas < fmaconf.gnumbases) error = malloc (nmeas * sizeof(complex float));

	if (!mpirank) fprintf (stderr, "Second DBIM pass (%d iterations)\n", dbimit[1]);
	for (dbimit[1] += i; i < dbimit[1]; ++i) {
		errnorm = cgnorm = 0;
		if (nmeas < fmaconf.gnumbases) {
			/* Find the minimum-norm solution. */
			errnorm = dbimerr (error, NULL, field, &hislv, &hislv, &srcmeas, &obsmeas);
			cgnorm = cgmn (error, crt, &hislv, &srcmeas, &obsmeas, regparm[1]);
		} else {
			/* Find the least-squares solution. */
			errnorm = dbimerr (NULL, rn, field, &hislv, &hislv, &srcmeas, &obsmeas);
			cgnorm = cgls (rn, crt, &hislv, &srcmeas, &obsmeas, regparm[1]);
		}
		
		for (j = 0; j < fmaconf.numbases; ++j)
			fmaconf.contrast[j] += crt[j];

		sprintf (fname, "%s.inverse.%03d", outproj, i);
		prtcontrast (fname, fmaconf.contrast,
				fmaconf.bslist, fmaconf.numbases, gsize);

		if (!mpirank)
			fprintf (stderr, "DBIM: %g, CG: %g (%d).\n", errnorm, cgnorm, i);
		if (errnorm < tolerance[1]) break;
	}

	ScaleME_finalizeParHostFMA ();

	freedircache ();

	free (fmaconf.contrast);
	free (fmaconf.radpats);
	free (field);
	if (error) free (error);
	free (srcmeas.locations);
	free (obsmeas.locations);
	delfrechbuf ();

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

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
#include "direct.h"
#include "mlfma.h"
#include "io.h"
#include "cg.h"

#include "util.h"

void usage (char *name) {
	fprintf (stderr, "Usage: %s [-s #] [-r #] [-a #] [-o <output prefix>] -i <input prefix>\n", name);
	fprintf (stderr, "\t-i <input prefix>: Specify input file prefix\n");
	fprintf (stderr, "\t-o <output prefix>: Specify output file prefix (defaults to input prefix)\n");
	fprintf (stderr, "\t-s #: The number of per-leap simultaneous views (default: 1)\n");
	fprintf (stderr, "\t-r #: The number of iterations for spectral radius estimation (default: none)\n");
	fprintf (stderr, "\t-a: Use ACA with specified tolerance for far-field transformations\n");
}

float dbimerr (complex float *error, complex float *rn, complex float *field,
	solveparm *hislv, solveparm *loslv, measdesc *src, measdesc *obs) {
	complex float *rhs, *crt, *err, *fldptr;
	float errnorm = 0, lerr, errd = 0;
	int j, i, nit;
	long k, nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;

	if (!error) err = malloc (obs->count * sizeof(complex float));
	else err = error;

	rhs = malloc (2L * nelt * sizeof(complex float));
	crt = rhs + nelt;

	if (rn) memset (rn, 0, nelt * sizeof(complex float));

	for (j = 0, fldptr = field; j < src->count; ++j, fldptr += obs->count) {
		/* Build the right-hand side for the specified location. Use
		 * point sources, rather than plane waves, for excitation. */
		buildrhs (rhs, src->locations + 3 * j);

		MPI_Barrier (MPI_COMM_WORLD);
		/* Run the iterative solver. The solution is stored in the RHS. */
		for (i = 0, nit = 1; i < hislv->restart && nit > 0; ++i)
			nit = bicgstab (rhs, crt, i, hislv->maxit, hislv->epscg, 1);

		/* Convert total field into contrast current. */
		for (k = 0; k < nelt; ++k)
			crt[k] *= fmaconf.contrast[k];

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

/* Establish the regularization parameter using an adapation of the
 * Lavarello and Oelze method described in their 2009 DBIM paper. */
float scalereg (float fact, float mse) {
	return fact * mse * mse * mse;
}

int main (int argc, char **argv) {
	char ch, *inproj = NULL, *outproj = NULL, **arglist, fname[1024];
	int mpirank, mpisize, i, nmeas, dbimit[2], q, stride = 1,
	    gsize[3], specit = 0, useaca = 0;
	complex float *rn, *crt, *field, *fldptr, *error, *refct;
	float errnorm = 0, tolerance[2], regparm[4], erninc,
	      trange[2], prange[2], crtmse = 0.0, gamma, sigma = 1.0;
	solveparm hislv, loslv;
	measdesc obsmeas, srcmeas, ssrc;
	long nelt, j;
	float acatol = -1;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpisize);

	if (!mpirank) fprintf (stderr, "Square-cell acoustic MLFMA.\n");

	arglist = argv;

	while ((ch = getopt (argc, argv, "i:o:s:r:a:")) != -1) {
		switch (ch) {
		case 'i':
			inproj = optarg;
			break;
		case 'o':
			outproj = optarg;
			break;
		case's':
			stride = atoi(optarg);
			break;
		case 'r':
			specit = atoi(optarg);
			break;
		case 'a':
			acatol = atof(optarg);
			useaca = 1;
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

	/* Read the basic configuration for only one observation shell. */
	sprintf (fname, "%s.input", inproj);
	getconfig (fname, &hislv, &loslv, &srcmeas, &obsmeas, 1);

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
	ScaleME_preconf (useaca);
	ScaleME_getListOfLocalBasis (&(fmaconf.numbases), &(fmaconf.bslist));

	nelt = (long)fmaconf.numbases * (long)fmaconf.bspboxvol;

	/* Allocate the RHS vector, residual vector and contrast. */
	fmaconf.contrast = malloc (3L * nelt * sizeof(complex float));
	rn = fmaconf.contrast + nelt;
	crt = rn + nelt;

	/* Store the grid size for writing of contrast values. */
	gsize[0] = fmaconf.nx; gsize[1] = fmaconf.ny; gsize[2] = fmaconf.nz;

	if (!mpirank) fprintf (stderr, "Reading local portion of contrast file.\n");
	/* Read the guess contrast for the local basis set. */
	sprintf (fname, "%s.guess", inproj);
	/* First try the group-ordered file, then try the old matrix format. */
	getctgrp (fmaconf.contrast, fname, gsize,
			fmaconf.bslist, fmaconf.numbases, fmaconf.bspbox);

	/* Read the reference contrast for the local basis set. */
	if (!mpirank)  fprintf (stderr, "Reading local portion of reference file.\n");
	refct = malloc (nelt * sizeof(complex float));
	sprintf (fname, "%s.reference", inproj);
	/* If the reference file is not found, set the reference to NULL to avoid
	 * checking the MSE at each step. Try both formats. */
	if (!getctgrp (refct, fname, gsize, fmaconf.bslist,
				fmaconf.numbases, fmaconf.bspbox)) {
		free (refct);
		refct = NULL;
	}


	/* Precalculate some values for the FMM and direct interactions. */
	fmmprecalc (acatol, useaca);
	i = dirprecalc ();
	if (!mpirank) fprintf (stderr, "Finished precomputing %d near interactions.\n", i);

	/* Finish the ScaleME initialization. */
	ScaleME_postconf ();

	/* Allocate the observation array. */
	field = malloc (nmeas * sizeof(complex float));

	/* Read the measurements and compute their norm. */
	getfields (inproj, field, obsmeas.count, srcmeas.count, &erninc);

	/* Build the root interpolation matrix for measurements. */
	ScaleME_buildRootInterpMat (obsmeas.imat, 6, obsmeas.ntheta,
			obsmeas.nphi, obsmeas.trange, obsmeas.prange);

	/* Build the root interpolation matrix for adjoint Frechet fields. Note the
	 * angular shift since the incoming fields have to be flipped. */
	trange[0] = M_PI - obsmeas.trange[0];
	trange[1] = M_PI - obsmeas.trange[1];
	prange[0] = M_PI + obsmeas.prange[0];
	prange[1] = M_PI + obsmeas.prange[1];
	ScaleME_buildRootInterpMat (obsmeas.imat + 1, 6,
			obsmeas.ntheta, obsmeas.nphi, trange, prange);

	MPI_Barrier (MPI_COMM_WORLD);

	if (!mpirank) fprintf (stderr, "Initialization complete.\n");

	/* Set the number of per-leap transmissions. */
	ssrc.count = stride = MAX(1,stride);
	/* Allocate the transmitter specifications. */
	ssrc.locations = malloc(3 * ssrc.count * sizeof(float));
	/* Blank the interpolation matrices. */
	ssrc.imat[0] = ssrc.imat[1] = NULL;

	/* The error storage vector for a single transmitter group. */
	error = malloc (ssrc.count * obsmeas.count * sizeof(complex float));

	/* Calculate the spectral radius, if desired. */
	if (specit > 0 && dbimit[0] > 0) {
		if (!mpirank) fprintf (stderr, "Spectral radius estimation (%d iterations)\n", specit);
		sigma = specrad (specit, &loslv, &srcmeas, &obsmeas);
		if (!mpirank) fprintf (stderr, "Spectral radius: %g\n", sigma);
		regparm[0] *= sigma;
		regparm[1] *= sigma;
		regparm[2] *= sigma;
	}

	/* Start a two-pass DBIM with low tolerances first. */
	if (!mpirank) fprintf (stderr, "First DBIM pass (%d iterations)\n", dbimit[0]);
	for (i = 0, gamma = regparm[0]; i < dbimit[0]; ++i) {
		/* Use the leapfrogging (Kaczmarz-like) method. */
		for (q = 0; q < srcmeas.count; q += stride) {
			fldptr = field + q * obsmeas.count;

			/* Don't use more transmitters than available. */
			ssrc.count = MIN(srcmeas.count - q, stride);

			/* Set the transmitter locations. */
			memcpy (ssrc.locations, srcmeas.locations + 3 * q, 3 * ssrc.count * sizeof(float));
			errnorm = dbimerr (error, NULL, fldptr, &hislv, &loslv, &ssrc, &obsmeas);

			/* Solve the system with CG for minimum norm. */
			cgmn (error, crt, &loslv, &ssrc, &obsmeas, gamma);

			/* Update the background. */
			for (j = 0; j < nelt; ++j)
				fmaconf.contrast[j] += crt[j];

			if (refct) crtmse = mse (fmaconf.contrast, refct, nelt, 1);

			if (!mpirank)
				fprintf (stderr, "Sub-iteration %d/%d: RRE: %g, MSE: %g\n", i, q, errnorm, crtmse);
		}

		if (!mpirank) fprintf (stderr, "Reassess DBIM error.\n");
		errnorm = dbimerr (NULL, NULL, field, &hislv, &loslv, &srcmeas, &obsmeas);

		sprintf (fname, "%s.inverse.%03d", outproj, i);
		prtctgrp (fname, fmaconf.contrast, gsize, fmaconf.bslist,
				fmaconf.numbases, fmaconf.bspbox);

		if (refct) crtmse = mse (fmaconf.contrast, refct, nelt, 1);

		if (!mpirank) fprintf (stderr, "Iteration %d: RRE: %g, MSE: %g\n", i, errnorm, crtmse);
		if (errnorm < tolerance[0]) break;

		/* Skip regularization adjustment unless the skip count has
		 * been met or the parameter is already small enough. */
		if ((i + 1) % (int)(regparm[3]) || gamma <= regparm[1]) continue;

		/* Perform scaling of the regularization parameter. */
		gamma = scalereg(regparm[0], errnorm);
		if (!mpirank)
			fprintf (stderr, "Regularization parameter: %g\n", gamma);
	}

	free (error);
	error = NULL;

	/* Conduct non-leapfrog iterations now. */
	if (i < dbimit[0]) ++i;

	/* If the number of measurements is less than the number of pixels,
	 * use the minimum-norm conjugate gradient. */
	if (nmeas < fmaconf.gnumbases) error = malloc (nmeas * sizeof(complex float));

	/* Recalculate the spectral radius, if desired. */
	if (specit > 0 && dbimit[1] > 0) {
		regparm[0] /= sigma;
		regparm[1] /= sigma;
		regparm[2] /= sigma;
		if (!mpirank) fprintf (stderr, "Spectral radius estimation (%d iterations)\n", specit);
		sigma = specrad (specit, &loslv, &srcmeas, &obsmeas);
		if (!mpirank) fprintf (stderr, "Spectral radius: %g\n", sigma);
		regparm[0] *= sigma;
		regparm[1] *= sigma;
		regparm[2] *= sigma;
	}

	/* Set the regularization parameter for the next set of iterations. */
	gamma = scalereg(regparm[2], errnorm);

	if (!mpirank) fprintf (stderr, "Second DBIM pass (%d iterations)\n", dbimit[1]);
	for (dbimit[1] += i; i < dbimit[1]; ++i) {
		errnorm = 0;
		if (nmeas < fmaconf.gnumbases) {
			/* Find the minimum-norm solution. */
			errnorm = dbimerr (error, NULL, field, &hislv, &loslv, &srcmeas, &obsmeas);
			cgmn (error, crt, &loslv, &srcmeas, &obsmeas, gamma);
		} else {
			/* Find the least-squares solution. */
			errnorm = dbimerr (NULL, rn, field, &hislv, &loslv, &srcmeas, &obsmeas);
			cgls (rn, crt, &loslv, &srcmeas, &obsmeas, gamma);
		}

		for (j = 0; j < nelt; ++j) fmaconf.contrast[j] += crt[j];

		sprintf (fname, "%s.inverse.%03d", outproj, i);
		prtctgrp (fname, fmaconf.contrast, gsize, fmaconf.bslist,
				fmaconf.numbases, fmaconf.bspbox);

		if (refct) crtmse = mse (fmaconf.contrast, refct, nelt, 1);

		if (!mpirank) fprintf (stderr, "Iteration %d: RRE: %g, MSE: %g\n", i, errnorm, crtmse);
		if (errnorm < tolerance[1]) break;

		/* Skip regularization adjustment unless the skip count has
		 * been met or the parameter is already small enough. */
		if ((i + 1) % (int)(regparm[3]) || gamma <= regparm[1]) continue;

		/* Perform scaling of the regularization parameter. */
		gamma = scalereg(regparm[2], errnorm);
		if (!mpirank)
			fprintf (stderr, "Regularization parameter: %g\n", gamma);
	}

	ScaleME_finalizeParHostFMA ();

	freedircache ();

	free (fmaconf.contrast);
	free (fmaconf.radpats);
	free (field);
	if (error) free (error);
	delmeas (&srcmeas);
	delmeas (&obsmeas);
	delmeas (&ssrc);

	if (refct) free (refct);

	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}

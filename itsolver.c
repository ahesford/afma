#include <complex.h>
#include <string.h>

#include <mpi.h>

/* These headers are provided by ScaleME. */
#include <Complex.h>
#include <ScaleME.h>
#include <gmres.h>

#include "mlfma.h"
#include "itsolver.h"

solveparm solver;

static int compcrt (Complex *dst, Complex *src) {
	int i;
	complex float val, buf;

	for (i = 0; i < fmaconf.numbases; ++i) {
		buf = src[i].re + I * src[i].im;
		val = buf * fmaconf.contrast[i];
		dst[i].re = creal (val);
		dst[i].im = cimag (val);
	}

	return fmaconf.numbases;
}

static int augcrt (Complex *dst, Complex *src) {
	int i;

	for (i = 0; i < fmaconf.numbases; ++i) {
		dst[i].re = src[i].re - dst[i].re;
		dst[i].im = src[i].im - dst[i].im;
	}

	return fmaconf.numbases;
}

int cgmres (complex float *rhs, complex float *sol, int silent) {
	int icntl[7], irc[5], lwork, info[3], i, myRank;
	float rinfo[2], cntl[5], ldot[2], gdot[2];
	Complex *zwork, *solbuf, *tx, *ty, *tz, lzdot;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	/* Allocate memory for work array. */
	lwork = solver.restart * solver.restart +
		solver.restart * (fmaconf.numbases + 5) + 5 * fmaconf.numbases + 1;
	zwork = calloc (lwork, sizeof(Complex));
	solbuf = malloc (fmaconf.numbases * sizeof(Complex));

	/* Initialize the parameters. */
	initcgmres_(icntl, cntl);

	/* Only the root process should print convergence history. */
	if (!myRank && !silent) icntl[2] = 6;
	else icntl[2] = 0;

	/* None of the processes should make noise about GMRES. */
	if (silent) icntl[0] = icntl[1] = icntl[2] = 0;

	/* Decide if a preconditioner should be used. */
	if (solver.precond) icntl[3] = 1;
	else icntl[3] = 0;

	icntl[4] = 0; /* Use MGS for orthogonalization. */
	icntl[5] = 1; /* Use an initial guess: the incident field. */
	icntl[6] = solver.maxit; /* Set the maximum interation count. */

	cntl[0] = solver.epscg;

	/* copy the initial guess: use the RHS. */
	memcpy (zwork, rhs, fmaconf.numbases * sizeof(complex float));
	/* Copy the unpreconditioned RHS. */
	memcpy (zwork + fmaconf.numbases, rhs, fmaconf.numbases * sizeof(complex float));

	do {
		drivecgmres_(&(fmaconf.gnumbases), &(fmaconf.numbases),
				&(solver.restart), &lwork, zwork, irc,
				icntl, cntl, info, rinfo);
		if (!(info[0]) && !(irc[0])) break;

		switch (irc[0]) {
		case EXIT: case RIGHT_PRECOND: break;
		case MATVEC:
			   /* fortran indices start from 1 */
			   tx = zwork+irc[1] - 1;
			   ty = zwork+irc[3] - 1;
			   compcrt (solbuf, tx);
			   ScaleME_applyParFMA(REGULAR, solbuf, ty);
			   augcrt (ty, tx);
			   break;
		case LEFT_PRECOND:
			   tx = zwork+irc[1] - 1;
			   ty = zwork+irc[3] - 1;

			   /* No preconditioner is desired */
			   if (!solver.precond) {
				   memcpy (ty, tx, fmaconf.numbases * sizeof(Complex));
				   break;
			   }

			   ScaleME_applyParBP (0, tx, ty);
			   break;
		case DOT_PROD:
			   ty = zwork + irc[2] - 1;
			   tx = zwork + irc[1] - 1;
			   tz = zwork + irc[3] - 1;
			   for (i = 0; i < irc[4]; ++i) {
				   /* compute the local dot product first */
				   lzdot = Cdotc(fmaconf.numbases, tx, ty);
				   /* now do a global reduce to get the final answer */
				   ldot[0] = lzdot.re;
				   ldot[1] = lzdot.im;
				   gdot[0] = gdot[1] = 0.0;
				   MPI_Allreduce(ldot, gdot, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
				   tz->re = gdot[0];
				   tz->im = gdot[1];
				   ++tz;
				   tx += fmaconf.numbases;
			   }
			   break;
		}
	} while (irc[0]);

	if (info[0] && !myRank)
		fprintf (stdout, "CGMRES: return value: %d\n", info[0]); 

	memcpy (sol, zwork, fmaconf.numbases * sizeof(complex float));
	
	if (!myRank && !silent)
		fprintf(stdout, "CGMRES: %d iterations, %.6E PBE, %.6E BE.\n", info[1], rinfo[0], rinfo[1]);

	free (zwork);
	free (solbuf);
	return info[0];
}

#include <mpi.h>

/* These headers are provided by ScaleME. */
#include <Complex.h>
#include <ScaleME.h>
#include <gmres.h>

#include "mlfma.h"
#include "itsolver.h"

solveparm solver;

int cgmres (Complex *rhs, Complex *sol) {
	int icntl[7], irc[5], lwork, info[3], i, myRank;
	float rinfo[2], cntl[5], ldot[2], gdot[2];
	Complex *zwork = NULL, *tx, *ty, *tz, lzdot;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	/* Allocate memory for work array. */
	lwork = solver.restart * solver.restart +
		solver.restart * (fmaconf.numbases + 5) + 5 * fmaconf.numbases + 1;
	zwork = calloc (lwork, sizeof(Complex));

	/* Initialize the parameters. */
	initcgmres_(icntl, cntl);

	/* Only the root process should print convergence history. */
	if (!myRank) icntl[2] = 6;
	else icntl[2] = 0;

	/* Decide if a preconditioner should be used. */
	if (solver.precond) icntl[3] = 1;
	else icntl[3] = 0;

	icntl[4] = 0; /* Use MGS for orthogonalization. */
	icntl[5] = 0; /* Don't use an initial guess. */
	icntl[6] = solver.maxit; /* Set the maximum interation count. */
	icntl[7] = 0; /* Don't use an explicit MVP to compute residual. */

	cntl[0] = solver.epscg;

	/* copy the input and solution vectors to zwork */
	Ccopy(fmaconf.numbases, sol, zwork);
	Ccopy(fmaconf.numbases, rhs, zwork + fmaconf.numbases);

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
			   ScaleME_applyParFMA(REGULAR, tx, ty);
			   break;
		case LEFT_PRECOND:
			   tx = zwork+irc[1] - 1;
			   ty = zwork+irc[3] - 1;
			   
			   /* No preconditioner is desired */
			   if (!solver.precond) {
				   Ccopy(fmaconf.numbases, tx, ty);
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

	if (info[0]) fprintf (stdout, "CGMRES: return value: %d\n", info[0]); 
	
	Ccopy(fmaconf.numbases, zwork, sol);
	
	if (!myRank) {
		fprintf(stdout, "CGMRES: %d iterations\n", info[1]);
		fprintf(stdout, "CGMRES: workspace size: %d\n", info[2]);
		fprintf(stdout, "CGMRES: Preconditioned B.E.: %.6E\n", rinfo[0]);
		fprintf(stdout, "CGMRES: Non-preconditioned B.E.: %.6E\n", rinfo[1]);
	}

	free(zwork);
	return info[0];
}

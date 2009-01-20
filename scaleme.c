/* 21-jan-2006 hesford; incorporate initialization routines into fastant */
/* 22-jun-2000 sanjay; added the shared memory model */
/* 14-dec-1999 sanjay; made the following modifications:
	1. removed ScaleME_setComLevel call.
	2. replaced the call to ScaleME_initSeqHostFMA with
	   ScaleME_initParHostDataStructs as a part of parallelizing the code
	   for parallel GMRES
	3. similar changes have been made to finalize calls.
	4. removed initialization of PDAX routines. */
/* 19-aug-1998 sanjay; created */
/* 03-aug-1998 sanjay; port C.C.Lu and Weng Chew's TRIMOM code to ScaleME */
/* scaleme.c : interface to scaleme */	

/*
scaleme.c:
Created by Sanjay Velamparambil on 03-August-1998.
Copyright: Sanjay Velamparambil, Weng Cho Chew, University of Illinois.
*/ 

#include <stdio.h>
#include <mpi.h>

#include <ScaleME.h> /* Provided by the ScaleME library. */
#include "scaleme.h" /* The local ScaleME include, defining these functions. */

#include "mlfma.h"
#include "itsolver.h"

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (void) {
	int error;

	/* The problem and tree are both three-dimensional. */
	ScaleME_setDimen (3);
	ScaleME_setTreeType (3);

	/* The fields are scalar-valued. */
	ScaleME_setFields (1);

	/* The wave number is real-valued. */
	ScaleME_setWaveNumber (fmaconf.k0);

	/* Set some MLFMA parameters. */
	ScaleME_setNumBasis (fmaconf.gnumbases);
	ScaleME_setMaxLevel (fmaconf.maxlev);
	ScaleME_setPrecision (fmaconf.precision);
	ScaleME_setMAC (fmaconf.numbuffer);
	ScaleME_setInterpOrder (fmaconf.interpord);

	ScaleME_sqSetPreSorted (1);
	ScaleME_setTopComputeLevel (fmaconf.toplev);
	ScaleME_useSharedMemoryModel (fmaconf.sharedmax);

	/* Disk storage is not allowed. */
	ScaleME_useDiskStorage (NULL, 0);
	
	if (fmaconf.smallbox > 0)
		ScaleME_setSmallestBoxSize(fmaconf.smallbox);

	/*  let all processes start the initialisation together */
	MPI_Barrier(MPI_COMM_WORLD); 
	
	/* open the std files */
	MPFMA_stdout = stdout;
	MPFMA_stderr = stderr; 
	
	error = ScaleME_initSetUp (MPI_COMM_WORLD, interaction,
			radpattern, rcvpattern, bscenter);

	if (error) {
		fprintf(stdout, "ERROR: ScaleME pre-init failed.\n");
		ScaleME_finalizeParHostFMA();
		return -1;
	}

	return 0;
}

int ScaleME_postconf (void) {
	long mem = 0;

	if (ScaleME_completeSetUp()) {
		fprintf(stdout, "ERROR: ScaleME setup routine failed.\n");
		goto error_handle;
	}

	if (ScaleME_initParHostDataStructs()) {
		fprintf(stdout, "ERROR: ScaleME parallel host init failed.\n");
		goto error_handle;
	}

	if (solver.precond) {
		if (ScaleME_initParBP (&mem)) {
			fprintf (stdout, "ERROR: ScaleME BD preconditioner failed.\n");
			goto error_handle;
		}
	}

	return 0;

error_handle:
	ScaleME_finalizeParHostFMA();
	return -1;
}

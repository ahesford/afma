#include <string.h>
#include <complex.h>

#include <mpi.h>

#include <Complex.h> // ScaleME-provided complex include. */
#include <ScaleME.h>

#include "measure.h"
#include "mlfma.h"
#include "fsgreen.h"

measdesc measures;

int farfield (complex float *currents, int nthetas, int nphis,
		float *thetas, float *phis, complex float *result) {
	Complex *fields = (Complex *)result;
	/* The result must be a pointer to the pointer. */
	if (ScaleME_evlRootFarFld_PI (6, nthetas, nphis, thetas, phis,
				(Complex *)currents, &fields)) return 0;

	return 1;
}

int directfield (complex float *currents, complex float *result) {
	int i, j, gi;
	complex float val, *buf;
	float *cen;

	buf = malloc (measures.numobs * sizeof(complex float));
	memset (buf, 0, measures.numobs * sizeof(complex float));

	/* Compute the fields radiated by the local fields. */
	for (j = 0; j < fmaconf.numbases; ++j) {
		gi = fmaconf.bslist[j];
		bscenter (gi, cen);
		for (i = 0; i < measures.numobs; ++i) {
			val = fsgreen (fmaconf.k0, cen, measures.obsloc + 3 * i);
			buf[i] += val * currents[j];
		}
	}

	/* Multiply by cell volume. */
	for (i = 0; i < measures.numobs; ++i)
		buf[i] *= fmaconf.cell[0] * fmaconf.cell[1] * fmaconf.cell[2];

	/* Collect the solution across all processors. */
	MPI_Allreduce (result, buf, 2 * measures.numobs, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return measures.numobs;
}

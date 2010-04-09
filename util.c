#include <complex.h>
#include <math.h>

#include <mpi.h>

#include "util.h"

float mse (complex float *test, complex float *ref, long n) {
	long i;
	float err[2] = {0.0, 0.0}, tmp;

	/* Calculate local contributions to the mean-squared error. */
	for (i = 0; i < n; ++i) {
		tmp = cabs(test[i] - ref[i]);
		err[0] += tmp * tmp;
		tmp = cabs(ref[i]);
		err[1] += tmp * tmp;
	}

	/* Sum the local contributions. */
	MPI_Allreduce (MPI_IN_PLACE, err, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return sqrt(err[0] / err[1]);
}

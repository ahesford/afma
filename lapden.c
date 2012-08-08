#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <complex.h>
#include <math.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "precision.h"
#include "util.h"

/* Print the usage. */
void usage (char *name) {
	fprintf (stderr, "USAGE: %s [-h] <-d d> [-r r] [input [output]]\n", name);
	fprintf (stderr, "\t-h: Display this message and exit\n");
	fprintf (stderr, "\t-d: Use a sample spacing d in wavelengths (default: 0.1)\n");
	fprintf (stderr, "\t-r: Use a reference density r (default: 1000)\n");
	fprintf (stderr, "\tInput file name may be '-' or omitted for stdin\n");
	fprintf (stderr, "\tOutput file name may be '-' or omitted for stdout\n");
}

/* Compute the i-th frequency bin for an m-point DFT with spacing h. */
real fftfreq(int i, int m, real h) {
	int half = (m - 1) / 2;
	real f;

	f = 2. * M_PI * (real) ((i <= half) ? i : (i - m)) / (h * (real) m);
	return f;
}

/* Compute the density contrast term. */
long lapden (real d, real r, FILE *input, FILE *output) {
	int rsize[3], csize[3];
	long p, pc, m;
	real *rdat;
	cplx *cdat;
	FFTW_PLAN fplan, bplan;

	/* Read the 3-D matrix size (nx, ny, nz) from the input. */
	fread (rsize, sizeof(int), 3, input);
	p = (long) rsize[0] * (long) rsize[1] * (long) rsize[2];

	/* The most-rapidly-varying dimension of the output is shorter. */
	csize[0] = rsize[0] / 2 + 1;
	csize[1] = rsize[1];
	csize[2] = rsize[2];
	pc = (long) csize[0] * (long) csize[1] * (long) csize[2];

	/* Allocate the storage arrays. */
	rdat = FFTW_MALLOC (p * sizeof(real));
	cdat = FFTW_MALLOC (pc * sizeof(cplx));

	fprintf (stderr, "INFO: Planning %d x %d x %d FFTs\n", rsize[0], rsize[1], rsize[2]);
	/* Plan the FFTs. The array dimensions must be transposed to agree
	 * with the FORTRAN order of the input and output files. */
	fplan = FFTW_PLAN_DFT_R2C_3D (rsize[2], rsize[1], rsize[0],
			rdat, cdat, FFTW_ESTIMATE);
	bplan = FFTW_PLAN_DFT_C2R_3D (rsize[2], rsize[1], rsize[0],
			cdat, rdat, FFTW_ESTIMATE);

	/* Read the 3-D matrix into the real data array. */
	fread (rdat, sizeof(real), p, input);

	/* Compute the reciprocal of the square root of the relative density. */
#pragma omp parallel for default(shared) private(m)
	for (m = 0; m < p; ++m) rdat[m] = sqrt(r / rdat[m]);

	/* Compute the DFT of the density function. */
	FFTW_EXECUTE(fplan);

#pragma omp parallel default(shared)
{
		long i, j, k, l;
		real kx, ky, kz;

		/* Now calculate the scaled Laplacian term. */
#pragma omp for
		for (l = 0; l < pc; ++l) {
			/* Pull out the three-dimensional index. */
			i = l % csize[0];
			j = (l / csize[0]) % csize[1];
			k = l / (csize[0] * csize[1]);

			/* Compute the DFT frequency bins for the index. */
			kx = fftfreq(i, rsize[0], d);
			ky = fftfreq(j, rsize[1], d);
			kz = fftfreq(k, rsize[2], d);

			/* Scale by the spectral Laplacian factor.
			 * Also scale by array dimensions to counter FFT scaling. */
			cdat[l] *= -(kx * kx + ky * ky + kz * kz) / (real) p;
		}
}

	/* Compute the inverse DFT of the density function. */
	FFTW_EXECUTE(bplan);

	/* Write the output to the output file. */
	fwrite (rsize, sizeof(int), 3, output);
	fwrite (rdat, sizeof(real), p, output);

	FFTW_FREE (rdat);
	FFTW_FREE (cdat);
	FFTW_DESTROY_PLAN (fplan);
	FFTW_DESTROY_PLAN (bplan);
	return p;
}

int main (int argc, char **argv) {
	char ch, *progname;
	FILE *input = NULL, *output = NULL;
	real d = 0.1, r = 1000.0;

	/* Store the name used to invoke the program. */
	progname = argv[0];

	/* Process the input arguments. */
	while ((ch = getopt (argc, argv, "hd:r:")) != -1) {
		switch (ch)  {
		case 'd':
			d = (real) strtod (optarg, NULL);
			break;
		case 'r':
			r = (real) strtod (optarg, NULL);
			break;
		default:
			usage (progname);
			exit (EXIT_FAILURE);
		}
	}

	/* Point argv to the input and output specifications. */
	argc -= optind;
	argv += optind;

	/* Use stdin or open an input file. */
	if (argc < 1 || !strcmp("-", argv[0])) input = stdin;
	else if (!(input = fopen(argv[0], "rb"))) {
		fprintf (stderr, "ERROR: Could not open %s\n", argv[0]);
		exit (EXIT_FAILURE);
	}

	/* Use stdout or open an output file. */
	if (argc < 2 || !strcmp("-", argv[1])) output = stdout;
	else if (!(output = fopen(argv[1], "wb"))) {
		fprintf (stderr, "ERROR: Could not open %s\n", argv[1]);
		exit (EXIT_FAILURE);
	}

#ifdef _OPENMP
	/* Configure FFTW to use threads, if possible. */
	int nt;
	nt = FFTW_INIT_THREADS();
	if (nt == 0) {
		fprintf (stderr, "ERROR: Could not initialize FFTW threads\n");
		exit (EXIT_FAILURE);
	}
	nt = omp_get_max_threads();
	FFTW_PLAN_WITH_NTHREADS(nt);
	fprintf (stderr, "INFO: FFTW will use %d threads\n", nt);
#endif

	lapden (d, r, input, output);

#ifdef _OPENMP
	FFTW_CLEANUP_THREADS();
#else
	FFTW_CLEANUP();
#endif

	fclose(input);
	fclose(output);

	return EXIT_SUCCESS;
}

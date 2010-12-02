#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>

#include "util.h"

/* Print the usage. */
void usage (char *name) {
	fprintf (stderr, "USAGE: %s <-n bpg | -r [-t Nx,Ny,Nz]> [input [output]]\n", name);
	fprintf (stderr, "\t-n: Number of basis functions per group per dimension\n");
	fprintf (stderr, "\t-r: Reverse the mapping to build a matrix file\n");
	fprintf (stderr, "\t-t: Truncate the matrix when reverse mapping is used\n");
	fprintf (stderr, "\tInput file name may be '-' or omitted for stdin\n");
	fprintf (stderr, "\tOutput file name may be '-' or omitted for stdout\n");
}

/* Convert an input matrix file into a grouped file
 * with bpg basis functions per group per dimension. */
int mat2grp (FILE *matfile, FILE *grpfile, int bpg) {
	int msize[3], grphdr[5], bpgvol, rrow;
	long i, grpslab, matslab;
	complex float *grpvals, *slabvals;

	bpgvol = bpg * bpg * bpg;

	/* Read the matrix size. */
	fread (msize, sizeof(int), 3, matfile);

	/* Set and write the group header. */
	grphdr[0] = 0;
	grphdr[1] = GEDIV(msize[0], bpg);
	grphdr[2] = GEDIV(msize[1], bpg);
	grphdr[3] = GEDIV(msize[2], bpg);
	grphdr[4] = bpg;
	fwrite (grphdr, sizeof(int), 5, grpfile);

	matslab = msize[0] * msize[1];
	grpslab = grphdr[1] * grphdr[2] * bpgvol;

	grpvals = malloc (grpslab * sizeof(complex float));
	slabvals = malloc (bpg * matslab * sizeof(complex float));

	for (i = 0; i < msize[2]; i += bpg) {
		/* Read a slab of data. Watch not to run off the end. */
		rrow = MIN(bpg, MAX(0, msize[2] - i));
		fread (slabvals, sizeof(complex float), rrow * matslab, matfile);

		/* Blank the slab of groups. */
		memset (grpvals, 0, grpslab * sizeof(complex float));

#pragma omp parallel default(shared)
{
		long l, j, k;
		int gidx[3], bidx[3];

		/* Now build each of the groups in the slab. */
#pragma omp for
		for (j = 0; j < grpslab; ++j) {
			/* Pull out the group index and the element index. */
			l = j / bpgvol;
			k = j % bpgvol;

			/* Find the index of the starting basis of the group. */
			GRID(gidx, l, grphdr[1], grphdr[2]);
			gidx[0] *= bpg;
			gidx[1] *= bpg;

			/* Find the index of the basis in the slab. */
			GRID(bidx, k, bpg, bpg);
			bidx[0] += gidx[0];
			bidx[1] += gidx[1];

			/* Skip values off the end of the matrix. */
			if (bidx[0] >= msize[0] || bidx[1] >= msize[1]
					|| bidx[2] >= rrow) continue;

			k = bidx[0] + msize[0] * (bidx[1] + msize[1] * bidx[2]);

			/* Copy the value into position. */
			grpvals[j] = slabvals[k];
		}
}

		/* Write the slab of groups. */
		fwrite (grpvals, sizeof(complex float), grpslab, grpfile);
	}

	free (grpvals);
	free (slabvals);
	return grphdr[1] * grphdr[2] * grphdr[3];
}

/* Convert an input grouped file into an output matrix file. */
int grp2mat (FILE *grpfile, FILE *matfile, int *mtrunc) {
	int msize[3], grphdr[5], bpgvol, rrow, bpg;
	long i, grpslab, matslab;
	complex float *grpvals, *slabvals;

	/* Read the group header:
	 * 0 Nx Ny Nz Bpg
	 * (Nx,Ny,Nz) is size of grid of groups.
	 * Bpg is number of basis functions per group per dimension. */
	fread (grphdr, sizeof(int), 5, grpfile);

	/* The number of elements in a group per dimension. */
	bpg = grphdr[4];
	/* The total number of elements in a group. */
	bpgvol = bpg * bpg * bpg;

	/* Build the overall matrix size. */
	msize[0] = grphdr[1] * bpg;
	msize[1] = grphdr[2] * bpg;
	msize[2] = grphdr[3] * bpg;

	/* Don't truncate the matrix if not desired. */
	if (!mtrunc) mtrunc = msize;

	/* The number of matrix elements in one z-slab. */
	matslab = mtrunc[0] * mtrunc[1];
	/* The number of matrix elements in one z-slab of groups. */
	grpslab = grphdr[1] * grphdr[2] * bpgvol;

	/* Allocate a slab of data for rearranging. */
	grpvals = malloc (grpslab * sizeof(complex float));
	slabvals = malloc (bpg * matslab * sizeof(complex float));

	/* Write the matrix size header. */
	fwrite (mtrunc, sizeof(int), 3, matfile);

	/* Loop through slabs of groups to do the conversion. */
	for (i = 0; i < mtrunc[2]; i += bpg) {
		/* Track the number of output rows to write for this slab. */
		rrow = MIN(bpg, MAX(0, mtrunc[2] - i));

		/* Read the slab of groups. */
		fread (grpvals, sizeof(complex float), grpslab, grpfile);

		/* Blank the matrix slab. */
		memset (slabvals, 0, rrow * matslab * sizeof(complex float));

#pragma omp parallel default(shared)
{
		long l, j, k;
		int gidx[3], bidx[3];

		/* Now build each slab. */
#pragma omp for
		for (j = 0; j < grpslab; ++j) {
			/* Pull out the group and element index. */
			l = j / bpgvol;
			k = j % bpgvol;

			/* Find the index of the starting basis of the group. */
			GRID(gidx, l, grphdr[1], grphdr[2]);
			gidx[0] *= bpg;
			gidx[1] *= bpg;

			/* Find the index of the basis in the slab. */
			GRID(bidx, k, bpg, bpg);
			bidx[0] += gidx[0];
			bidx[1] += gidx[1];

			/* Don't write values that should be truncated. */
			if (bidx[0] >= mtrunc[0] || bidx[1] >= mtrunc[1]
					|| bidx[2] >= rrow) continue;

			k = bidx[0] + mtrunc[0] * (bidx[1] + mtrunc[1] * bidx[2]);

			/* Copy the value into position. */
			slabvals[k] = grpvals[j];
		}
}

		/* Write the slab of the matrix. */
		fwrite (slabvals, sizeof(complex float), rrow * matslab, matfile);
	}

	free (grpvals);
	free (slabvals);
	return mtrunc[0] * mtrunc[1] * mtrunc[2];
}

int main (int argc, char **argv) {
	char ch, *progname;
	int bpg = 0, rev = 0, trunc[3], usetrunc = 0;
	FILE *input = NULL, *output = NULL;

	/* Store the name used to invoke the program. */
	progname = argv[0];

	/* Process the input arguments. */
	while ((ch = getopt (argc, argv, "rn:t:h")) != -1) {
		switch (ch)  {
		case 'r':
			rev = 1;
			break;
		case 'n':
			bpg = strtol(optarg, NULL, 0);
			break;
		case 't':
			usetrunc = 1;
			trunc[0] = strtol(strtok(optarg, ","), NULL, 0);
			trunc[1] = strtol(strtok(NULL, ","), NULL, 0);
			trunc[2] = strtol(strtok(NULL, ","), NULL, 0);
		default:
			usage (progname);
			exit (EXIT_FAILURE);
		}
	}

	/* Point argv to the input and output specifications. */
	argc -= optind;
	argv += optind;

	/* Respect the mutual exclusivity of -n and -r, and require one.
	 * Also ensure that an input and output file have both been specified. */
	if ((rev && bpg > 0) || (!rev && bpg < 1)) {
		usage (progname);
		exit (EXIT_FAILURE);
	}

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

	if (!rev) mat2grp (input, output, bpg);
	else grp2mat (input, output, usetrunc ? trunc : NULL);

	fclose(input);
	fclose(output);

	return EXIT_SUCCESS;
}

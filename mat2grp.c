#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>

#include "util.h"

/* Print the usage. */
void usage (char *name) {
	printf ("USAGE: %s <-n bpg | -r [-t Nx,Ny,Nz]> <-i input> <-o output>\n", name);
	printf ("\t-n: Number of basis functions per group per dimension\n");
	printf ("\t-r: Reverse the mapping to build a matrix file\n");
	printf ("\t-t: Truncate the matrix when reverse mapping is used\n");
	printf ("\t-i: Input file name (default: stdin)\n");
	printf ("\t-o: Output file name (default: stdout)\n");
}

/* Convert an input matrix file into a grouped file
 * with bpg basis functions per group per dimension. */
int mat2grp (FILE *matfile, FILE *grpfile, int bpg) {
	int msize[3], grphdr[5], gidx[3], bidx[3], k, bpgvol, rrow;
	long offset, i, ntg;
	complex float *grpvals;

	bpgvol = bpg * bpg * bpg;

	grpvals = malloc (bpgvol * sizeof(complex float));

	/* Read the matrix size. */
	fread (msize, sizeof(int), 3, matfile);

	/* Set and write the group header. */
	grphdr[0] = 0;
	grphdr[1] = GEDIV(msize[0], bpg);
	grphdr[2] = GEDIV(msize[1], bpg);
	grphdr[3] = GEDIV(msize[2], bpg);
	grphdr[4] = bpg;
	fwrite (grphdr, sizeof(int), 5, grpfile);

	/* The total number of groups. */
	ntg = (long)grphdr[1] * (long)grphdr[2] * (long)grphdr[3];

	/* Now loop through all groups, reading the values
	 * from matfile and writing to the group file. */
	for (i = 0; i < ntg; ++i) {
		/* Get the global, 3-D group index. */
		GRID(gidx, i, grphdr[1], grphdr[2]);
		/* Find the starting basis index from the group index. */
		gidx[0] *= bpg;
		gidx[1] *= bpg;
		gidx[2] *= bpg;

		/* Don't try to read beyond the bounds of the file. */
		rrow = MIN(bpg, MAX(0, msize[0] - gidx[0]));

		/* Zero the group buffer. */
		memset (grpvals, 0, bpgvol * sizeof(complex float));

		for (k = 0; k < bpgvol; k += bpg) {
			/* Compute the index of the basis within the group. */
			GRID(bidx, k, bpg, bpg);
			bidx[0] += gidx[0];
			bidx[1] += gidx[1];
			bidx[2] += gidx[2];

			/* Don't read this row if it exceeds the grid bounds. */
			if (bidx[1] >= msize[1] || bidx[2] >= msize[2]) continue;

			/* Compute the position of the row in the file. */
			offset = 3L * sizeof(int) + sizeof(complex float) *
				((long)bidx[0] + (long)msize[0] * (long)bidx[1] +
				 (long)msize[0] * (long)msize[1] * (long)bidx[2]);
			fseek (matfile, offset, SEEK_SET);
			fread (grpvals + k, sizeof(complex float), rrow, matfile);
		}

		/* Write the group values in one block. */
		fwrite (grpvals, sizeof(complex float), bpgvol, grpfile);
	}

	free (grpvals);
	return ntg;
}

/* Convert an input grouped file into an output matrix file. */
int grp2mat (FILE *grpfile, FILE *matfile, int *mtrunc) {
	int msize[3], grphdr[5], gidx[3], bidx[3], k, bpgvol, rrow, fd;
	long i, ntg, offset;
	complex float *grpvals;

	/* Read the group header:
	 * 0 Nx Ny Nz Bpg
	 * (Nx,Ny,Nz) is size of grid of groups.
	 * Bpg is number of basis functions per group per dimension. */
	fread (grphdr, sizeof(int), 5, grpfile);

	/* Note the total number of groups. */
	ntg = (long)grphdr[1] * (long)grphdr[2] * (long)grphdr[3];

	/* Allocate a read buffer for each group. */
	bpgvol = grphdr[4] * grphdr[4] * grphdr[4];
	grpvals = malloc (bpgvol * sizeof(complex float));

	/* Build the overall matrix size. */
	msize[0] = grphdr[1] * grphdr[4];
	msize[1] = grphdr[2] * grphdr[4];
	msize[2] = grphdr[3] * grphdr[4];

	/* Don't truncate the matrix if not desired. */
	if (!mtrunc) mtrunc = msize;

	/* Truncate the output file and write the header. */
	fd = fileno(matfile);
	offset = 3L * sizeof(int) + ntg * (long)bpgvol * sizeof(complex float);
	ftruncate (fd, offset);
	fwrite (mtrunc, sizeof(int), 3, matfile);

	/* Now loop through all groups, reading the values
	 * from matfile and writing to the group file. */
	for (i = 0; i < ntg; ++i) {
		/* Get the global, 3-D group index. */
		GRID(gidx, i, grphdr[1], grphdr[2]);
		/* Find the starting basis index from the group index. */
		gidx[0] *= grphdr[4];
		gidx[1] *= grphdr[4];
		gidx[2] *= grphdr[4];

		/* Read the group values in one block. */
		fread (grpvals, sizeof(complex float), bpgvol, grpfile);

		/* Don't try to read beyond the bounds of the file. */
		rrow = MIN(grphdr[4], MAX(0, mtrunc[0] - gidx[0]));

		/* Write the values of the group. */
		for (k = 0; k < bpgvol; k += grphdr[4]) {
			/* Compute the index of the basis within the group. */
			GRID(bidx, k, grphdr[4], grphdr[4]);
			bidx[0] += gidx[0];
			bidx[1] += gidx[1];
			bidx[2] += gidx[2];

			/* Don't read this row if it exceeds the grid bounds. */
			if (bidx[1] >= mtrunc[1] || bidx[2] >= mtrunc[2]) continue;

			/* Compute the position of the row in the file. */
			offset = 3L * sizeof(int) + sizeof(complex float) *
				((long)bidx[0] + (long)mtrunc[0] * (long)bidx[1] +
				 (long)mtrunc[0] * (long)mtrunc[1] * (long)bidx[2]);
			fseek (matfile, offset, SEEK_SET);
			fwrite (grpvals + k, sizeof(complex float), rrow, matfile);
		}
	}

	free (grpvals);
	return ntg;
}

int main (int argc, char **argv) {
	char ch, *inname = NULL, *outname = NULL;
	int bpg = 0, rev = 0, trunc[3], usetrunc = 0;
	FILE *input = NULL, *output = NULL;

	/* Process the input arguments. */
	while ((ch = getopt (argc, argv, "rn:i:o:t:")) != -1) {
		switch (ch)  {
		case 'r':
			rev = 1;
			break;
		case 'n':
			bpg = strtol(optarg, NULL, 0);
			break;
		case 'i':
			inname = optarg;
			break;
		case 'o':
			outname = optarg;
			break;
		case 't':
			usetrunc = 1;
			trunc[0] = strtol(strtok(optarg, ","), NULL, 0);
			trunc[1] = strtol(strtok(NULL, ","), NULL, 0);
			trunc[2] = strtol(strtok(NULL, ","), NULL, 0);
		default:
			usage (argv[0]);
			exit (EXIT_FAILURE);
		}
	}

	/* Respect the mutual exclusivity of -n and -r, and require one.
	 * Also ensure that an input and output file have both been specified. */
	if ((rev && bpg > 0) || (!rev && bpg < 1) || !inname || !outname) {
		usage (argv[0]);
		exit (EXIT_FAILURE);
	}

	if (!(input = fopen(inname, "rb"))) {
		fprintf (stderr, "ERROR: Could not open %s\n", inname);
		exit (EXIT_FAILURE);
	}

	if (!(output = fopen(outname, "wb"))) {
		fprintf (stderr, "ERROR: Could not open %s\n", outname);
		exit (EXIT_FAILURE);
	}

	if (!rev) mat2grp (input, output, bpg);
	else grp2mat (input, output, usetrunc ? trunc : NULL);

	fclose(input);
	fclose(output);

	return EXIT_SUCCESS;
}

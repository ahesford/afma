#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <complex.h>
#include <float.h>
#include <math.h>

#define ALVTOL FLT_EPSILON
#define PTEQ(a,b) ((a)[0] == (b)[0] && (a)[1] == (b)[1] && (a)[2] == (b)[2])

typedef struct {
	float c, a, r;
} material;

/* Compute the wave number for dimensionless distances given the sound speed
 * (m/s), attenuation (dB / cm / MHz), and background sound speed (m/s). */
complex float wavenum (float c, float a, float cb) {
	float kr, ki;

	/* Unitless imaginary part. */
	ki = 1e-4 * cb * a * log(10) / 20;
	/* Relative real part. */
	kr = 2 * M_PI * cb / c;

	return kr + I * ki;
}

/* Compute the Laplacian of the density at a specified point. */
float lapden (float *r, float *lr, float *nr, float delta,
		long pos, long nx, long ny) {
	long x, y;
	float dlap, nval, pval;

	x = pos % nx;
	y = pos / nx;

	/* The self contribution to the Laplacian. */
	dlap = -6.0 / r[pos];

	/* Contribution of the x offsets with bounds checking. */
	if (x >= nx - 1) nval = 1.0;
	else nval = 1.0 / r[pos + 1];

	if (x <= 0) pval = 1.0;
	else pval = 1.0 / r[pos - 1];

	dlap += pval + nval;

	/* Contribution of the y offsets with bounds checking. */
	if (y >= ny - 1) nval = 1.0;
	else nval = 1.0 / r[pos + nx];

	if (y <= 0) pval = 1.0;
	else pval = 1.0 / r[pos - nx];

	dlap += pval + nval;

	/* Contribution of the z offsets with bounds checking. */
	if (!nr) nval = 1.0;
	else nval = 1.0 / nr[pos];

	if (!lr) pval = 1.0;
	else pval = 1.0 / lr[pos];

	dlap += pval + nval;

	/* Divide out the step size. */
	dlap /= (delta * delta);

	return dlap;
}

/* Build the contrast of the current slab. Requires the previous and next
 * slab densities to evaluate the Laplacian of the density. */
long ctslab (complex float *contrast, complex float *k, float *r, float *lr,
		float *nr, complex float kbg, long nx, long ny, float delta) {
	long i, npx;
	float rval;

	npx = nx * ny;

#pragma omp parallel for default(shared) private(rval,i)
	for (i = 0; i < npx; ++i) {
		/* Compute the relative wave number. */
		contrast[i] = k[i] / kbg;
		/* Now square and shift to form the contrast. */
		contrast[i] = contrast[i] * contrast[i] - 1;

		/* Compute the Laplacian of the density. */
		rval = lapden (r, lr, nr, delta, i, nx, ny);

		/* Add in the density term. */
		contrast[i] -= rval * r[i] / (kbg * kbg);
	}

	return npx;
}


/* Build the next z-slab of the model for specified parameters. */
long modelslab (complex float *k, float *r, FILE *template, FILE *avlrand,
		FILE *adprand, long npx, material *bgval, material *ctval,
		material *avlval, material *adpval) {
	long i;
	float *tmap, *arnd, *frnd, cs = 1.0, as = 0.0, rs = 1.0;

	/* Allocate arrays to store the slab values. */
	tmap = malloc (3 * npx * sizeof(float));
	arnd = tmap + npx;
	frnd = arnd + npx;

	/* Read the tissue parameters for the current slab. */
	fread (tmap, sizeof(float), npx, template);
	fread (arnd, sizeof(float), npx, avlrand);
	fread (frnd, sizeof(float), npx, adprand);

#pragma omp parallel for default(shared) private(i,cs,as,rs)
	for (i = 0; i < npx; ++i) {
		switch((int)(tmap[i])) {
		/* Connective tissue values should be used here. */
		case 0:
			/* The dimensionless wave number. */
			k[i] = wavenum(ctval->c, ctval->a, bgval->c);
			/* Square root of the relative density. */
			r[i] = sqrt(ctval->r / bgval->r);
			break;
		/* Background values should be used here. */
		case 1:
			/* The dimensionless wave number. */
			k[i] = wavenum(bgval->c, bgval->a, bgval->c);
			/* Relative density in the background is always unity. */
			r[i] = 1.0;
			break;
		/* Random models should be used here. */
		default:
			if (arnd[i] > ALVTOL) {
				/* Use the alveolar models. */
				cs = avlval[0].c + avlval[1].c * arnd[i];
				as = avlval[0].a + avlval[1].a * arnd[i];
				rs = avlval[0].r + avlval[1].r * arnd[i];
			} else {
				/* Use the adipose models. */
				cs = adpval[0].c + adpval[1].c * frnd[i];
				as = adpval[0].a + adpval[1].a * frnd[i];
				rs = adpval[0].r + adpval[1].r * frnd[i];
			}

			/* The dimensionless wave number. */
			k[i] = wavenum(cs, as, bgval->c);
			/* Square root of the relative density. */
			r[i] = sqrt(rs / bgval->r);
		}
	}

	free (tmap);
	return npx;
}

/* Build the entire model, writing slab-by-slab to the output file. */
long buildmodel (FILE *output, FILE *template, FILE *avlrand, FILE *adprand,
		float delta, material *bgval, material *ctval,
		material *avlval, material *adpval) {
	complex float *kslab, *contrast, *k, *nk, kbg;
	float *rstore, *r, *lr, *nr;
	long npx, i, nx, ny, nz;
	int size[3], nsize[3];

	/* Build the backgroudn wave number. */
	kbg = wavenum(bgval->c, bgval->a, bgval->c);

	/* Read the size of each of the model grids. */
	fread (size, sizeof(int), 3, template);
	fread (nsize, sizeof(int), 3, avlrand);
	if (!PTEQ(nsize,size)) {
		fprintf (stderr, "ERROR: alveolar and template sizes do not agree.\n");
		return 0;
	}
	fread (nsize, sizeof(int), 3, adprand);
	if (!PTEQ(nsize,size)) {
		fprintf (stderr, "ERROR: adipose and template sizes do not agree.\n");
		return 0;
	}

	/* Write the size of the model grid to the contrast output. */
	fwrite (size, sizeof(int), 3, output);

	nx = size[0]; ny = size[1]; nz = size[2];
	npx = nx * ny;

	fprintf (stderr, "Grid size is %ld x %ld x %ld.\n", nx, ny, nz);

	/* Allocate the arrays used for storage. */
	rstore = malloc (3 * npx * sizeof(float));
	kslab = malloc (3 * npx * sizeof(complex float));
	contrast = kslab + 2 * npx;

	/* Point to the data stores. */
	k = kslab;
	nk = kslab + npx;

	lr = NULL;
	r = rstore;
	nr = rstore + npx;

	/* Build the model slab for the first plane. */
	modelslab (k, r, template, avlrand, adprand, npx, bgval,
			ctval, avlval, adpval);

	for (i = 1; i < nz; ++i) {
		if (!(i % 10)) fprintf (stderr, "Building model for slab %ld\n", i);

		/* Read the next slab of the material. */
		modelslab (nk, nr, template, avlrand, adprand, npx, bgval,
				ctval, avlval, adpval);

		/* Build and write the contrast for the previously-read slab. */
		ctslab (contrast, k, r, lr, nr, kbg, nx, ny, delta);
		fwrite (contrast, sizeof(complex float), npx, output);

		/* Update the media pointers. */
		k = kslab + (i % 2) * npx;
		nk = kslab + ((i + 1) % 2) * npx;

		lr = rstore + ((i - 1) % 3) * npx;
		r = rstore + (i % 3) * npx;
		nr = rstore + ((i + 1) % 3) * npx;
	}

	/* Build and write out the last slab. */
	ctslab (contrast, k, r, lr, NULL, kbg, nx, ny, delta);
	fwrite (contrast, sizeof(complex float), npx, output);

	free (kslab);
	free (rstore);

	return npx;
}

/* Print the usage. */
void usage (char *name) {
	printf ("USAGE: %s [-b <c,a,r>] [-c <c,a,r>] -f <c,a,r,c,a,r> -a <c,a,r,c,a,r>\n", name);
	printf ("\t\t-d <step> -t <template> -r <alveolar> -R <adipose> -o <contrast>\n");
	printf ("\t-b: Background sound speed, attenuation, and density (default: water)\n");
	printf ("\t-c: Connective tissue sound speed, attenuation, and density (default: 1545,0,1120)\n");
	printf ("\t-f: Sound speed, attenuation, and density (and swings)\n");
	printf ("\t    for adipose (fat) random tisse model\n");
	printf ("\t-a: Sound speed, attenuation, and density (and swings)\n");
	printf ("\t    for alveolar random tisse model\n");
	printf ("\t-d: Grid step size in wavelengths (default: 0.1)\n");
	printf ("\t-t: Tissue template file name\n");
	printf ("\t-r: Alveolar random file name\n");
	printf ("\t-R: Adipose random file name\n");
	printf ("\t-o: Output contrast file name\n");
}

/* Parse the string token defining tissue types. */
int parsetype (material *mat, int nmat, char *str) {
	char *val, *sep = ",";
	int i;

	/* Initialize to the first token. */
	if (!(val = strtok (str, sep))) return -1;

	/* Keep breaking up the string to assign values to the materials. */
	for (i = 0; i < nmat; ++i) {
		mat[i].c = atof (val);
		if (!(val = strtok (NULL, sep))) return -1;
		mat[i].a = atof (val);
		if (!(val = strtok (NULL, sep))) return -1;
		mat[i].r = atof (val);
		if (!(val = strtok (NULL, sep))) return -1;
	}

	return 0;
}

int main (int argc, char **argv) {
	char **arglist, ch, *tmplname = NULL, *avlname = NULL,
	     *adpname = NULL, *outname = NULL;
	material bgval, ctval, avlval[2], adpval[2];
	float delta = 0.1;
	FILE *template, *avlrand, *adprand, *output;

	arglist = argv;

	/* Set defaults for the background and connective tissue. */
	bgval.c = 1509; bgval.a = 0; bgval.r =  997;
	ctval.c = 1545; ctval.a = 0; ctval.r = 1120;

	/* Process the input arguments. */
	while ((ch = getopt (argc, argv, "b:c:f:a:d:t:r:R:o:h")) != -1) {
		switch (ch)  {
		case 'b':
			parsetype (&bgval, 1, optarg);
			break;
		case 'c':
			parsetype (&ctval, 1, optarg);
			break;
		case 'f':
			parsetype (adpval, 2, optarg);
			break;
		case 'a':
			parsetype (avlval, 2, optarg);
			break;
		case 'd':
			delta = atof (optarg);
			break;
		case 't':
			tmplname = optarg;
			break;
		case 'r':
			avlname = optarg;
			break;
		case 'R':
			adpname = optarg;
			break;
		case 'o':
			outname = optarg;
			break;
		default:
			usage (arglist[0]);
			exit (EXIT_FAILURE);
		}
	}

	/* Make sure all file arguments were provided. */
	if (tmplname == NULL || avlname == NULL || adpname == NULL || outname == NULL) {
		usage (arglist[0]);
		exit (EXIT_FAILURE);
	}

	/* Open all files. */
	if (!(template = fopen (tmplname, "r"))) {
		fprintf (stderr, "ERROR: could not open template file %s\n", tmplname);
		exit (EXIT_FAILURE);
	}
	if (!(avlrand = fopen (avlname, "r"))) {
		fprintf (stderr, "ERROR: could not open alveolar file %s\n", avlname);
		exit (EXIT_FAILURE);
	}
	if (!(adprand = fopen (adpname, "r"))) {
		fprintf (stderr, "ERROR: could not open adipose file %s\n", adpname);
		exit (EXIT_FAILURE);
	}
	if (!(output = fopen (outname, "w"))) {
		fprintf (stderr, "ERROR: could not open output file %s\n", outname);
		exit (EXIT_FAILURE);
	}

	/* Build the model slice by slice. */
	buildmodel (output, template, avlrand, adprand, delta,
			&bgval, &ctval, avlval, adpval);

	/* Close all open files. */
	fclose (output);
	fclose (adprand);
	fclose (avlrand);
	fclose (template);

	return EXIT_SUCCESS;
}

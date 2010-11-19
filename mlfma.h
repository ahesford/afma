#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include "util.h"

typedef struct {
	float min[3], cen[3], cell, cellvol, grplen;
	float precision;
	int nx, ny, nz, gnumbases, numbases;
	int bspbox, maxlev, numbuffer, interpord, toplev, bspboxvol;
	int fo2iterm, fo2iord, fo2iosr;
	int *bslist, nsamp, acarank;
	float k0;
	complex float *contrast, *radpats;
} fmadesc;

extern fmadesc fmaconf;

static inline void bscenter (int gi, float *cen) {
	int idx[3];

	GRID(idx, gi, fmaconf.nx, fmaconf.ny);

	cen[0] = fmaconf.min[0] + ((float)idx[0] + 0.5) * fmaconf.grplen;
	cen[1] = fmaconf.min[1] + ((float)idx[1] + 0.5) * fmaconf.grplen;
	cen[2] = fmaconf.min[2] + ((float)idx[2] + 0.5) * fmaconf.grplen;
}

void acafarpattern (int, int *, void *, void *, float *, int);
void farpattern (int, int *, void *, void *, float *, int);

int fmmprecalc (float, int);

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (int);
int ScaleME_postconf (void);

#endif /* __MLFMA_H_ */

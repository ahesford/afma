#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h>
#include <fftw3.h>

#include "mlfma.h"

typedef struct {
	float min[3], cen[3], cell, cellvol;
	float precision;
	int nx, ny, nz, gnumbases, numbases;
	int bspbox, maxlev, numbuffer, interpord, toplev;
	int fo2iterm, fo2iord, fo2iosr;
	int *bslist;
	float k0;
	complex float *contrast, *radpats, *gridints;
	int nsamp, nbors, nborsex;
	fftwf_plan fplan, bplan;
} fmadesc;

extern fmadesc fmaconf;

static inline void bsindex (int gi, int *idx) {
	idx[0] = gi % fmaconf.nx;
	idx[1] = (gi / fmaconf.nx) % fmaconf.ny;
	idx[2] = gi / (fmaconf.nx * fmaconf.ny);
}

static inline void bscenter (int gi, float *cen) {
	int idx[3];

	bsindex (gi, idx);

	cen[0] = fmaconf.min[0] + ((float)idx[0] + 0.5) * fmaconf.cell;
	cen[1] = fmaconf.min[1] + ((float)idx[1] + 0.5) * fmaconf.cell;
	cen[2] = fmaconf.min[2] + ((float)idx[2] + 0.5) * fmaconf.cell;
}

void farpattern (int, int *, void *, void *, float *, int);
int fmmprecalc ();
void blockinteract (int, int, int *, int *, void *, void *, float *);

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (void);
int ScaleME_postconf (void);

#endif /* __MLFMA_H_ */

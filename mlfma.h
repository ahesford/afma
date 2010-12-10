#ifndef __MLFMA_H_
#define __MLFMA_H_

#include "precision.h"
#include "util.h"

typedef struct {
	real min[3], cen[3], cell, cellvol, grplen;
	real precision;
	int nx, ny, nz, gnumbases, numbases;
	int bspbox, maxlev, numbuffer, interpord, toplev, bspboxvol;
	int fo2itxlev, fo2ibclev, fo2iord, fo2iosr;
	int *bslist, nsamp, acarank;
	real k0;
	cplx *contrast, *radpats;
} fmadesc;

extern fmadesc fmaconf;

static inline void bscenter (int gi, real *cen) {
	int idx[3];

	GRID(idx, gi, fmaconf.nx, fmaconf.ny);

	cen[0] = fmaconf.min[0] + ((real)idx[0] + 0.5) * fmaconf.grplen;
	cen[1] = fmaconf.min[1] + ((real)idx[1] + 0.5) * fmaconf.grplen;
	cen[2] = fmaconf.min[2] + ((real)idx[2] + 0.5) * fmaconf.grplen;
}

void acafarpattern (int, int *, void *, void *, real *, int);
void farpattern (int, int *, void *, void *, real *, int);

int fmmprecalc (real, int);

/* initialisation and finalisation routines for ScaleME */
int ScaleME_preconf (int);
int ScaleME_postconf (void);

#endif /* __MLFMA_H_ */

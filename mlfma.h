#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h>
#include "mlfma.h"

typedef struct {
	float min[3], max[3], cell[3];
	float precision, smallbox;
	int nx, ny, nz, gnumbases, numbases;
	int maxlev, numbuffer, interpord, toplev, sharedmax;
	int *bslist;
	float k0;
	complex float *contrast;
	float nbors[3];
	complex float *gridints;
} fmadesc;

extern fmadesc fmaconf;

void radpattern (int, float *, float *, void *);
void rcvpattern (int, float *, float *, void *);

void impedance (int, int, void *);

static inline void bsindex (int gi, int *idx) {
	idx[0] = gi % fmaconf.nx;
	idx[1] = (gi / fmaconf.nx) % fmaconf.ny;
	idx[2] = gi / (fmaconf.nx * fmaconf.ny);
}

static inline void bscenter (int gi, float *cen) {
	int idx[3];

	bsindex (gi, idx);

	cen[0] = fmaconf.min[0] + ((float)idx[0] + 0.5) * fmaconf.cell[0];
	cen[1] = fmaconf.min[1] + ((float)idx[1] + 0.5) * fmaconf.cell[1];
	cen[2] = fmaconf.min[2] + ((float)idx[2] + 0.5) * fmaconf.cell[2];
}

void interaction (int, int, void *);

int preimpedance ();
void blockinteract (int, int, int *, int *, void *, void *);

#endif /* __MLFMA_H_ */

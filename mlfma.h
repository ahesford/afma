#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h>

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
void bscenter (int, float *);

void bsindex (int, int *);

void interaction (int, int, void *);

int preimpedance ();
void blockinteract (int, int, int *, int *, void *, void *);

#endif /* __MLFMA_H_ */

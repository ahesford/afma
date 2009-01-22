#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h> // This is provided with GCC.
#include <Complex.h> // This is provided by ScaleME.

typedef struct {
	float min[3], max[3], cell[3];
	float precision, smallbox;
	int nx, ny, nz, gnumbases, numbases;
	int maxlev, numbuffer, interpord, toplev, sharedmax;
	int *bslist;
	float k0;
	complex float *contrast;
} fmadesc;

extern fmadesc fmaconf;

void radpattern (int, float *, float *, Complex *);
void rcvpattern (int, float *, float *, Complex *);

void impedance (int, int, Complex *);
void bscenter (int, float *);

#endif /* __MLFMA_H_ */

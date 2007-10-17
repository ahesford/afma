#ifndef __MLFMA_H_
#define __MLFMA_H_

#include <complex.h> // This is provided with GCC.
#include <Complex.h> // This is provided by ScaleME.

typedef struct {
	float min[3], max[3], cell[3];
	int nx, ny, nz;
	complex float k0;
} fmadesc;

extern fmadesc mlfma;

void radpattern (int, float *, float *, Complex *);
void rcvpattern (int, float *, float *, Complex *);

void impedance (int, int, Complex *);
void bscenter (int, float *);

#endif /* __MLFMA_H_ */

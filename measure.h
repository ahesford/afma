#ifndef __MEASURE_H_
#define __MEASURE_H_

#include <Complex.h>

typedef struct {
	int numsrc, numobs;
	float *srcloc, *obsloc;
} measdesc;

extern measdesc measures;

int farfield (complex float *, int, int, float *, float *, complex float *);
int directfield (complex float *, complex float *);

#endif /* __MEASURE_H_ */

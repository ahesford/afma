#ifndef __MEASURE_H_
#define __MEASURE_H_

#include <complex.h>

typedef struct {
	int count, ntheta, nphi;
	float radius;
	float prange[2], trange[2], *locations;
	void *imat[2];
} measdesc;

int buildrhs (complex float *, float *);

int farfield (complex float *, measdesc *, complex float *);
int buildlocs (measdesc *);
void delmeas (measdesc *);

#endif /* __MEASURE_H_ */

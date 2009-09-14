#ifndef __MEASURE_H_
#define __MEASURE_H_

#include <ScaleME_Complex.h>

typedef struct {
	int count, ntheta, nphi;
	float radius;
	float prange[2], trange[2], *locations;
} measdesc;

int farfield (complex float *, measdesc *, complex float *);
int directfield (complex float *, measdesc *, complex float *, complex float *);
int buildlocs (measdesc *);

#endif /* __MEASURE_H_ */

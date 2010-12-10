#ifndef __MEASURE_H_
#define __MEASURE_H_

#include "precision.h"

typedef struct {
	int count, ntheta, nphi;
	real prange[2], trange[2], *locations;
	void *imat[2];
} measdesc;

int buildrhs (cplx *, real *);

int farfield (cplx *, measdesc *, cplx *);
int buildlocs (measdesc *);
void delmeas (measdesc *);

#endif /* __MEASURE_H_ */

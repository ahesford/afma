#ifndef __MEASURE_H_
#define __MEASURE_H_

#include "precision.h"

typedef struct {
	int count, ntheta, nphi, plane;
	real prange[2], trange[2], *locations;
	void *imat[2];
} measdesc;

int buildrhs (cplx *, real *, int, real *);

int farfield (cplx *, measdesc *, cplx *);

int buildsrc (measdesc *, char *);
int buildobs (measdesc *, char *);
void delmeas (measdesc *);

#endif /* __MEASURE_H_ */

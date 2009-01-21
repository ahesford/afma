#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include <complex.h>

complex float radint (float, float *, float *, float *, float *, int, double *, double *);
complex float srcint (float, float *, float *, float *, int, double *, double *);
complex float oneptint (float, float *, float *, float *);
complex float selfint (float, float *);
complex float selfinthigh (float, float *);

#endif /* __INTEGRATE_H_ */

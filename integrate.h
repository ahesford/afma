#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include <complex.h>

#define OUTWT 0.347854845137454
#define INWT 0.652145154862546

#define OUTPT 0.861136311594053
#define INPT 0.339981043584856

#define NUMPTS 4

complex float srcint (float, float *, float *, float *);
complex float selfint (float, float *);

#endif /* __INTEGRATE_H_ */

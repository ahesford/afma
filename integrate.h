#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include <complex.h>

typedef complex float (*ifunc)(float, float *, float *);

complex float rcvint (float, float *, float *, float, ifunc);
complex float srcint (float, float *, float *, float, ifunc);
complex float selfint (float, float, int);

void bldintrules (int, int);
void delintrules ();

#endif /* __INTEGRATE_H_ */

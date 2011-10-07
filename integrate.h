#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include "precision.h"

typedef cplx (*ifunc)(real, real *, real *);

cplx rcvint (real, real *, real *, real, ifunc);
cplx srcint (real, real *, real *, real, ifunc);
cplx duffyint (real, real);
cplx selfint (real, real);

void bldintrules (int, int);
void delintrules ();

#endif /* __INTEGRATE_H_ */

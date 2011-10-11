#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include "precision.h"

/* A function to serve as an integrand. */
typedef cplx (*ifunc)(real, real *, real *);
/* A function to serve as an integrator. */
typedef cplx (*integrator)(real, real *, real *, real *, ifunc);

/* A routine for double integration. */
cplx rcvint (real, real *, real *, real *, integrator, ifunc);

/* Routines for singular and smooth single integration. */
cplx srcint (real, real *, real *, real *, ifunc);
cplx duffyint (real, real *, real *, real *, ifunc);

void bldintrules (int, int);
void delintrules ();

#endif /* __INTEGRATE_H_ */

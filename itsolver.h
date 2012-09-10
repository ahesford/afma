#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include "precision.h"

typedef struct {
  int restart, maxit;
  real epscg;
} solveparm;

typedef struct {
	cplx *z, *az;
	int nmax, ntot, start;
} augspace;

int matvec (cplx *, cplx *, cplx *, int);
int gmres (cplx *, cplx *, int, int, real, int, augspace *);
int bicgstab (cplx *, cplx *, int, int, real, int);

#endif /* __ITSOLVER_H_ */

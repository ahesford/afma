#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <complex.h>

typedef struct {
  int restart, maxit;
  float epscg;
} solveparm;

typedef struct {
	complex float *z, *az;
	int nmax, ntot, start;
} augspace;

int matvec (complex float *, complex float *, complex float *);
int gmres (complex float *, complex float *, int, int, float, int, augspace *);
int bicgstab (complex float *, complex float *, int, int, float, int);

#endif /* __ITSOLVER_H_ */

#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <complex.h>

typedef struct {
  int restart, maxit;
  float epscg;
} solveparm;

int matvec (complex float *, complex float *, complex float *);
int cgmres (complex float *, complex float *, int, int, solveparm *);
int bicgstab (complex float *, complex float *, int, int, float, int);

#endif /* __ITSOLVER_H_ */

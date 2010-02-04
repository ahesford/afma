#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <complex.h>

typedef struct {
  int restart, maxit;
  float epscg;
} solveparm;

int cgmres (complex float *, complex float *, int, solveparm *);
int bicgstab (complex float *, complex float *, int, solveparm *);

#endif /* __ITSOLVER_H_ */

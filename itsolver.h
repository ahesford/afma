#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <complex.h>

typedef struct {
  int restart, maxit;
  float epscg;
} solveparm;

float cgmres (complex float *, complex float *, int, solveparm *);
float bicgstab (complex float *, complex float *, int, solveparm *);

#endif /* __ITSOLVER_H_ */

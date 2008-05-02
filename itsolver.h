#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <Complex.h> /* ScaleME complex include */

typedef struct {
  int restart, maxit, precond;
  float epscg;
} solveparm;

extern solveparm solver;

int cgmres (complex float *, complex float *, int);

#endif /* __ITSOLVER_H_ */

#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <Complex.h> /* ScaleME complex include */

typedef struct {
  int restart, maxit, precond;
  float epscg;
} solveparm;

extern solveparm solver;

int cgmres (Complex *, Complex *);

#endif /* __ITSOLVER_H_ */
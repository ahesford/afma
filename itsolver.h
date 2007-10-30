#ifndef __ITSOLVER_H_
#define __ITSOLVER_H_

#include <Complex.h> /* ScaleME complex include */

typedef struct {
  int restart, maxit, precond;
  float epscg, regparm;
} solveparm;

extern solveparm solver;

int cgmres (complex float *, complex float *);

#endif /* __ITSOLVER_H_ */

#ifndef __FRECHET_H_
#define __FRECHET_H_

#include "precision.h"

#include "measure.h"
#include "itsolver.h"

int frechet (cplx *, cplx *, cplx *, measdesc *, solveparm *);
int frechadj (cplx *, cplx *, cplx *, measdesc *, solveparm *);
real specrad (int, solveparm *, measdesc *, measdesc *);

#endif /* __FRECHET_H_ */

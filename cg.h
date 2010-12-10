#ifndef __CG_H_
#define __CG_H_

#include "precision.h"

#include "measure.h"
#include "itsolver.h"

real cgls (cplx *, cplx *, solveparm *, measdesc *, measdesc *, real);
real cgmn (cplx *, cplx *, solveparm *, measdesc *, measdesc *, real);

#endif /* __CG_H_ */

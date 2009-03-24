#ifndef __CG_H_
#define __CG_H_

#include <complex.h>

#include "measure.h"
#include "itsolver.h"

float cgls (complex float *, complex float *,
		solveparm *, measdesc *, measdesc *, float);
float cgmn (complex float *, complex float *,
		solveparm *, measdesc *, measdesc *, float);

#endif /* __CG_H_ */

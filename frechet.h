#ifndef __FRECHET_H_
#define __FRECHET_H_

#include <complex.h>

#include "measure.h"
#include "itsolver.h"

int frechet (complex float *, complex float *,
		complex float *, measdesc *, solveparm *);
int frechadj (complex float *, complex float *,
		complex float *, measdesc *, solveparm *);

#endif /* __FRECHET_H_ */

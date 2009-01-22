#ifndef __EXCITE_H_
#define __EXCITE_H_

#include <complex.h>
#include "measure.h"

complex float planerhs (int, float *);
complex float pointrhs (int, float *);
int buildrhs (complex float *, float *, int);
int multirhs (complex float *, measdesc *, complex float *, int);

#endif /* __EXCITE_H_ */

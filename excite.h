#ifndef __EXCITE_H_
#define __EXCITE_H_

#include <complex.h>
#include "measure.h"

int buildrhs (complex float *, float *, int);
int multirhs (complex float *, measdesc *, complex float *, int);

int precompgrf (measdesc *, complex float *, int);
int precomprhs (complex float *, measdesc *, complex float *, complex float *);

#endif /* __EXCITE_H_ */

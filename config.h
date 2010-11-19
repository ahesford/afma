#ifndef __IO_H_
#define __IO_H_

#include <stdio.h>

#include "itsolver.h"
#include "measure.h"

void skipcomments (FILE *);
void getdbimcfg (char *, int *, float *, float *);
void getconfig (char *, solveparm *, solveparm *, measdesc *, measdesc *, int);

#endif /* __IO_H_ */

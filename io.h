#ifndef __IO_H_
#define __IO_H_

#include <stdio.h>
#include <complex.h>

#include "itsolver.h"
#include "measure.h"

void skipcomments (FILE *);
void getdbimcfg (char *, int *, float *, float *);
void getconfig (char *, solveparm *, solveparm *, measdesc *, measdesc *);
int getcontrast (complex float *, char *, int *, int);
int prtcontrast (char *, complex float *, int *, int *, int);
int writefld (char *, measdesc *, complex float *);
int readfld (complex float *, char *, int);
int getfields (char *, complex float *, int, int, float *);

#endif /* __IO_H_ */

#ifndef __IO_H_
#define __IO_H_

#include <stdio.h>
#include <complex.h>
#include "measure.h"

void skipcomments (FILE *);
void getconfig (char *);
void getcontrast (char *, int *, int);

int prtcontrast (char *, complex float *);
int prtfield (char *, measdesc *, complex float *);
int getfield (char *, complex float *, int);

#endif /* __IO_H_ */

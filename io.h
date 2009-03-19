#ifndef __IO_H_
#define __IO_H_

#include <stdio.h>
#include <complex.h>
#include "measure.h"

void skipcomments (FILE *);
void getdbimcfg (char *, int *, float *, float *);
void getconfig (char *);
void getcontrast (char *, int *, int);

int prtcontrast (char *, complex float *);
int prtfldhdr (char *, measdesc *, measdesc *);
int appendfld (char *, measdesc *, complex float *);
int getfields (char *, complex float *, int, float *);

#endif /* __IO_H_ */

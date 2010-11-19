#ifndef __IO_H_
#define __IO_H_

#include <complex.h>

int getctgrp (complex float *, char *, int *, int *, int, int);
int prtctgrp (char *, complex float *, int *, int *, int, int);
int writefld (char *, int, int, complex float *);
int readfld (complex float *, char *, int);
int getfields (char *, complex float *, int, int, float *);

#endif /* __IO_H_ */

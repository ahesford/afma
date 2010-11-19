#ifndef __IO_H_
#define __IO_H_

#include <complex.h>

int getcontrast (complex float *, char *, int *, int *, int, int);
int prtcontrast (char *, complex float *, int *, int *, int, int);
int writefld (char *, int, int, complex float *);
int readfld (complex float *, char *, int);
int getfields (char *, complex float *, int, int, float *);

#endif /* __IO_H_ */

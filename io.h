#ifndef __IO_H_
#define __IO_H_

#include "precision.h"

int getctgrp (cplx *, char *, int *, int *, int, int);
int chkctprt (char *);
int prtctgrp (char *, cplx *, int *, int *, int, int);
int writefld (char *, int, int, cplx *);
int readfld (cplx *, char *, int);
int getfields (char *, cplx *, int, int, real *);

#endif /* __IO_H_ */

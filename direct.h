#ifndef __DIRECT_H_
#define __DIRECT_H_

#include "precision.h"

int greengrid (cplx *, int, int, real, real, int *);
int dirprecalc (int);

cplx *cacheboxrhs (int, int);
void blockinteract(int, int, int *, int *, int);

int mkdircache ();
void clrdircache ();
void freedircache ();

#endif /* __DIRECT_H_ */

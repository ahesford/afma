#ifndef __DIRECT_H_
#define __DIRECT_H_

#include <complex.h>

int greengrid (complex float *, int, int, float, float, int *, int, int);
int dirprecalc (int, int);

complex float *cacheboxrhs (int, int);
void blockinteract(int, int, int *, int *, int);

int mkdircache ();
void clrdircache ();
void freedircache ();

#endif /* __DIRECT_H_ */

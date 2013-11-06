#ifndef __UTIL_H_
#define __UTIL_H_

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define GRID(a,gi,nx,ny) { (a)[0] = (gi) % (nx); \
	(a)[1] = ((gi) / (nx)) % (ny); \
	(a)[2] = (gi) / ((nx) * (ny)); }
#define IDX(nc,i,j,k) ((i) + (nc) * ((j) + (nc) * (k)))

#define GEDIV(a,b) ((a) / (b) + ((a) % (b) == 0 ? 0 : 1))

/* Condition for selective reorthogonalization, suggested by
 * Giraud and Langou in CERFACS Tech. Report No. TR/PA/02/52,
 * "Robust selective Gram-Schmidt reorthogonalization". */
#define IMGS_L 0.99

#include "precision.h"

#ifdef _ATLAS
#define log2(a) (log(a) / log(2))
double complex cexp (double complex);
#endif /* _ATLAS */

real sinc (real);
real mse (cplx *, cplx *, long, int);

int sampcoords (real *, int, int, int);
int cellcoords (real *, int, int, real);

int inset (int, int *, int);
int maxind (cplx *, int, int *, int);

int cmgs (cplx *, cplx *, cplx *, long, int);
cplx pardot (cplx *, cplx *, long);
real parnorm (cplx *, long);

int gaussleg (real *, real *, int);
#endif /* __UTIL_H_ */

#ifndef __UTIL_H_
#define __UTIL_H_

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
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

#ifdef _FREEBSD
#define log2(a) (log(a) / log(2))
double complex cexp (double complex);
#endif /* _FREEBSD */

float sinc (float);
float mse (complex float *, complex float *, long, int);

int sampcoords (float *, int, float *, int, int);
int cellcoords (float *, int, int, float);

int inset (int, int *, int);
int maxind (complex float *, int, int *, int);

int cmgs (complex float *, complex float *, complex float *, long, int);
complex float pardot (complex float *, complex float *, long);
float parnorm (complex float *, long);
#endif /* __UTIL_H_ */

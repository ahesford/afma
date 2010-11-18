#ifndef __UTIL_H_
#define __UTIL_H_

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define IDX(nc,i,j,k) ((i) + (nc) * ((j) + (nc) * (k)))
#define GRID(a,nc,i) { (a)[0] = (i) % (nc); \
	(a)[1] = ((i) / (nc)) % (nc); \
	(a)[2] = (i) / ((nc) * (nc)); }

#define GEDIV(a,b) ((a) / (b) + ((a) % (b) == 0 ? 0 : 1))

/* Condition for selective reorthogonalization, suggested by
/* Giraud and Langou in CERFACS Tech. Report No. TR/PA/02/52,
 * "Robust selective Gram-Schmidt reorthogonalization". */
#define IMGS_L 0.99

/* Define some no-op OpenMP functions if OpenMP is not used. */
#ifndef _OPENMP
int omp_get_max_threads () { return 1; }
int omp_get_thread_num () { return 0; }

int omp_set_lock (omp_lock_t *x) { return 0; }
int omp_init_lock (omp_lock_t *x) { return 0; }
int omp_unset_lock (omp_lock_t *x) { return 0; }
#endif /* _OPENMP */

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

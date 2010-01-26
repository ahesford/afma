#ifndef __UTIL_H_
#define __UTIL_H_

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define IDX(ny,nz,i,j,k) ((k) + (nz) * ((j) + (ny) * (i)))
#define SQIDX(ny,i,j,k) IDX(ny,ny,i,j,k)

/* Define some no-op OpenMP functions if OpenMP is not used. */
#ifndef _OPENMP
int omp_get_max_threads () { return 1; }
int omp_get_thread_num () { return 0; }
#endif /* _OPENMP */

#ifdef _FREEBSD
#define log2(a) (log(a) / log(2))
#endif /* _FREEBSD */

#endif /* __UTIL_H_ */

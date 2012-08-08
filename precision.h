#ifndef __PRECISION_H_
#define __PRECISION_H_

#include <complex.h>
#include <float.h>
#include <math.h>

/* The file ScaleME_real.h defines the "real" data type.
 * The file ScaleME_Complex.h defines the "cplx" data type. */
#include "ScaleME_Complex.h"
#include "ScaleME_real.h"

/* Define a generic complex and real type for single or double precision. */

#ifdef DOUBLEPREC /* Use double precision. */
#define REAL_EPSILON DBL_EPSILON

#define FFTW_INIT_THREADS fftw_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define FFTW_CLEANUP_THREADS fftw_cleanup_threads
#define FFTW_CLEANUP fftw_cleanup
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define FFTW_FREE fftw_free
#define FFTW_MALLOC fftw_malloc
#define FFTW_EXECUTE_DFT fftw_execute_dft
#define FFTW_EXECUTE fftw_execute
#define FFTW_PLAN_DFT_3D fftw_plan_dft_3d
#define FFTW_PLAN_DFT_R2C_3D fftw_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftw_plan_dft_c2r_3d
#define FFTW_PLAN fftw_plan

#define TRSV cblas_ztrsv
#define GEMV cblas_zgemv
#define GEMM cblas_zgemm
#define DOTC_SUB cblas_zdotc_sub
#define DOTU_SUB cblas_zdotu_sub

#define ROT zrot_
#define LARTG zlartg_
#define GEQRF zgeqrf_
#define UNGQR zungqr_
#define GESVD zgesvd_

#else /* Use single precision. */
#define REAL_EPSILON FLT_EPSILON

#define FFTW_INIT_THREADS fftwf_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define FFTW_CLEANUP_THREADS fftwf_cleanup_threads
#define FFTW_CLEANUP fftwf_cleanup
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define FFTW_FREE fftwf_free
#define FFTW_MALLOC fftwf_malloc
#define FFTW_EXECUTE_DFT fftwf_execute_dft
#define FFTW_EXECUTE fftwf_execute
#define FFTW_PLAN_DFT_3D fftwf_plan_dft_3d
#define FFTW_PLAN_DFT_R2C_3D fftwf_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftwf_plan_dft_c2r_3d
#define FFTW_PLAN fftwf_plan

#define TRSV cblas_ctrsv
#define GEMV cblas_cgemv
#define GEMM cblas_cgemm
#define DOTC_SUB cblas_cdotc_sub
#define DOTU_SUB cblas_cdotu_sub

#define ROT crot_
#define LARTG clartg_
#define GEQRF cgeqrf_
#define UNGQR cungqr_
#define GESVD cgesvd_

#endif

#endif /* __PRECISION_H_ */

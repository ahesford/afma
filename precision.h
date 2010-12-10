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

#define FFTW_FREE fftw_free
#define FFTW_MALLOC fftw_malloc
#define FFTW_EXECUTE_DFT fftw_execute_dft
#define FFTW_PLAN_DFT_3D fftw_plan_dft_3d
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

#define FFTW_FREE fftwf_free
#define FFTW_MALLOC fftwf_malloc
#define FFTW_EXECUTE_DFT fftwf_execute_dft
#define FFTW_PLAN_DFT_3D fftwf_plan_dft_3d
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

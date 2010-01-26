#ifndef __DIRECT_H_
#define __DIRECT_H_

#include <complex.h>
#include <fftw3.h>

typedef struct {
	int nbox[3], nboxprod, boxmin[3], boxmax[3], totelts;
	int nbors, nborsvol, nfftprod, nfft[3];
	int *boxfill;
	complex float *boxrhs, *gridints;
	fftwf_plan fplan, bplan;
} dirdesc;

extern dirdesc dircache;

int greengrid (complex float *, int, int, float, float, int *);
int dirprecalc ();

void blockinteract(int, int, int *, int *, void *, void *, float *);

#endif /* __DIRECT_H_ */

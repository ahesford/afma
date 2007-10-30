#ifndef __FRECHET_H_
#define __FRECHET_H_

#include <complex.h>

complex float *bldfrechbuf (int);
void delfrechbuf (void);

int frechet (complex float *, complex float *, complex float *);
int frechadj (complex float *, complex float *, complex float *);

#endif /* __FRECHET_H_ */

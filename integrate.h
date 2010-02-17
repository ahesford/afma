#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_

#include <complex.h>

/* Values for three-point Gauss-Legendre quadrature. */
#define OTWT3 0.555555555555556
#define CNWT3 0.888888888888889

#define OTPT3 0.774596669241483

/* Values for four-point Gauss-Legendre quadrature. */
#define OTWT4 0.347854845137454
#define INWT4 0.652145154862546

#define OTPT4 0.861136311594053
#define INPT4 0.339981043584856

/* Values for five-point Gauss-Legendre quadrature. */
#define OTWT5 0.236926885056189
#define INWT5 0.478628670499366
#define CNWT5 0.568888888888889

#define OTPT5 0.906179845938664
#define INPT5 0.538469310105683

typedef complex float (*ifunc)(float, float *, float *);

complex float rcvint (float, float *, float *, float, ifunc);
complex float srcint (float, float *, float *, float, ifunc);
complex float selfint (float, float);

#endif /* __INTEGRATE_H_ */

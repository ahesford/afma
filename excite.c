#include <math.h>
#include <complex.h>

#include "mlfma.h"
#include "integrate.h"
#include "fsgreen.h"
#include "excite.h"

/* Computes the RHS for a plane-wave in a given direction. */
complex float planerhs (int gi, float *srcdir) {
	complex float ans;
	float ctr[3];

	/* Find the center of the requested basis. */
	bscenter (gi, ctr);
	ans = rcvint (fsplane, fmaconf.k0, srcdir, ctr, fmaconf.cell);

	return ans;
}

/* Computes the RHS for a point source at the given location. */
complex float pointrhs (int gi, float *srcloc) {
	complex float ans;
	float ctr[3];

	bscenter (gi, ctr);
	ans = rcvint (fsgreen, fmaconf.k0, srcloc, ctr, fmaconf.cell);

	return ans;
}

/* Computes the entire RHS vector. */
int buildrhs (complex float *rhs, float *srcloc, int type) {
	int i;
	complex float (*rhsfunc) (int, float *);

	if (type) rhsfunc = planerhs;
	else rhsfunc = pointrhs;

	for (i = 0; i < fmaconf.numbases; ++i)
		rhs[i] = rhsfunc (fmaconf.bslist[i], srcloc);

	return fmaconf.numbases;
}

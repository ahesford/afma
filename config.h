#ifndef __CONFIG_H_
#define __CONFIG_H_

#include <stdio.h>

#include "precision.h"

#include "itsolver.h"

void skipcomments (FILE *);
void getdbimcfg (char *, int *, real *, real *);
void getconfig (char *, solveparm *, solveparm *);

#endif /* __CONFIG_H_ */

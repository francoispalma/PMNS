#ifndef PMNS4096_H
#define PMNS4096_H

#include "structs.h"

void pmns4096_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

#ifndef PMNS512_H
#define PMNS512_H

#include "structs.h"

void pmns512_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

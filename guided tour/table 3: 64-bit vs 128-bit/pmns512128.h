#ifndef PMNS512128_H
#define PMNS512128_H

#include "structs.h"

void pmns512128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);

#endif

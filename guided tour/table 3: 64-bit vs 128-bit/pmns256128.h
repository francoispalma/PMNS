#ifndef PMNS256128_H
#define PMNS256128_H

#include "structs.h"

void pmns256128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);

#endif

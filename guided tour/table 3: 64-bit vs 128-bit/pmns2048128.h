#ifndef PMNS2048128_H
#define PMNS2048128_H

#include "structs.h"

void pmns2048128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);

#endif

#include "pmns512128.h"
#include "params512128.h"
#include "corepmns128.c"

void pmns512128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// 512-bit modular multiplication using the polynomial reduction algorithm.
	pmns128_montg_mult(res, A, B);
}


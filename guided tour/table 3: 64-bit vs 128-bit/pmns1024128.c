#include "pmns1024128.h"
#include "params1024128.h"
#include "corepmns128.c"

void pmns1024128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// 1024-bit modular multiplication using the polynomial reduction algorithm.
	pmns128_montg_mult(res, A, B);
}


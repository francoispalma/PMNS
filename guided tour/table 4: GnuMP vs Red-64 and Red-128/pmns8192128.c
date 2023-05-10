#include "pmns8192128.h"
#include "params8192128.h"
#include "corepmns128.c"

void pmns8192128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// 8192-bit modular multiplication using the polynomial reduction algorithm.
	pmns128_montg_mult(res, A, B);
}


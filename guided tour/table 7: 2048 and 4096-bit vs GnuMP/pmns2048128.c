#include "pmns2048128.h"
#include "params2048128.h"
#include "corepmns128.c"

void pmns2048128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// 2048-bit modular multiplication using the polynomial reduction algorithm.
	pmns128_montg_mult(res, A, B);
}

void pmns2048128_montg_mult_2core(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// 2048-bit modular multiplication using the polynomial reduction algorithm.
	ppmns128_montg_mult(res, A, B, 2);
}

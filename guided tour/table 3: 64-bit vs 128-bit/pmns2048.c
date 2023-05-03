#include "pmns2048.h"
#include "params2048.h"
#include "corepmns.c"

void pmns2048_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 2048-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}


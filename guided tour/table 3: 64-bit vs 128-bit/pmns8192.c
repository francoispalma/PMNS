#include "pmns8192.h"
#include "params8192.h"
#include "corepmns.c"

void pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 8192-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}


#include "pmns1024.h"
#include "params1024.h"
#include "corepmns.c"

void pmns1024_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 1024-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}


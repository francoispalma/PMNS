#include "pmns512.h"
#include "params512.h"
#include "corepmns.c"

void pmns512_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 512-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}

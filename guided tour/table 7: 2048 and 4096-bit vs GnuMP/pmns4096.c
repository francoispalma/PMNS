#include "pmns4096.h"
#include "params4096.h"
#include "corepmns.c"

void pmns4096_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 4096-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}

void pmns4096_montg_mult_6core(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 4096-bit modular multiplication using the polynomial reduction algorithm.
	ppmns_montg_mult(res, A, B, 6);
}

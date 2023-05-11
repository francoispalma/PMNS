#include "pmns4096.h"
#include "params4096.h"
#include "corepmns.c"

void poly_pmns4096_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 4096-bit modular multiplication using the polynomial reduction algorithm.
	poly_pmns_montg_mult(res, A, B);
}

void latt_pmns4096_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 4096-bit modular multiplication using the lattice reduction algorithm.
	latt_pmns_montg_mult(res, A, B);
}


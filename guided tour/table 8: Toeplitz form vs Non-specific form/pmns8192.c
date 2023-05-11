#include "pmns8192.h"
#include "params8192.h"
#include "corepmns.c"

void poly_pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 8192-bit modular multiplication using the polynomial reduction algorithm.
	poly_pmns_montg_mult(res, A, B);
}

void latt_pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 8192-bit modular multiplication using the lattice reduction algorithm.
	latt_pmns_montg_mult(res, A, B);
}


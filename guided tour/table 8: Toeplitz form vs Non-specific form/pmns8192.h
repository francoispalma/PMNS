#ifndef PMNS8192_H
#define PMNS8192_H

#include "structs.h"

void poly_pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
void latt_pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

#ifndef PMNS256_H
#define PMNS256_H

#include "structs.h"

void poly_pmns256_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
void latt_pmns256_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

#ifndef PMNS1024_H
#define PMNS1024_H

#include "structs.h"

void poly_pmns1024_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
void latt_pmns1024_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

#ifndef PMNS8192_H
#define PMNS8192_H

#include "structs.h"

void pmns8192_montg_mult_5core(restrict poly res, const restrict poly A,
	const restrict poly B);
void pmns8192_montg_mult_6core(restrict poly res, const restrict poly A,
	const restrict poly B);
void pmns8192_montg_mult_8core(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

#ifndef PMNS2048_H
#define PMNS2048_H

#include "structs.h"

void pmns2048_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
void pmns2048_montg_mult_2core(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

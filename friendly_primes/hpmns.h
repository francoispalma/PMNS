#ifndef HPMNS_H_INCLUDED
#define HPMNS_H_INCLUDED

#include "structs.h"

void pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

void convert_binary_to_pmns(restrict poly res, const uint64_t* restrict op);

#endif

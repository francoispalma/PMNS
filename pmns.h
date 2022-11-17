#ifndef MONTGOM_H
#define MONTGOM_H

#include <stdint.h>
#include <stdlib.h>

#include "structs.h"
#include "params.h"

void convert_string_to_amns(restrict poly res, const char* string);
void convert_amns_to_multipre(restrict mpnum* res, const restrict poly P);

extern void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B);
extern void m1_mns_mod_mult_ext_red(int64_t* restrict R,
	const restrict poly A);
extern void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
extern void UNROLLED_amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
extern void mns_montg_int_red(restrict poly res, __int128* R);
extern void amns_sqandmult(restrict poly res, const char* restrict base,
	const char* restrict exponent);

#endif

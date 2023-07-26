#ifndef MPPMNS_H
#define MPPMNS_H

#include <stdint.h>
#include <stdlib.h>

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

#include "structs.h"
#include "params128.h"


void convert_string_to_amns128(restrict poly128 res, const char* string);
void convert_amns128_to_multipre(restrict mpnum* res, const restrict poly128 P);

void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B);
void m1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A);
void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
void UNROLLED_amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
void amns128_sqandmult(restrict poly128 res, const char* restrict base,
	const char* restrict exponent);

#endif

#ifndef HPMNS_H
#define HPMNS_H

#include <stdint.h>
#include "structs.h"

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B);
/*void m_pmns_mod_mult_ext_red(__int128* restrict R,
	const int64_t* restrict A);
void m1_pmns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A);*/
void pmns_montg_int_red(restrict poly res, __int128* restrict R);
void pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
void sea512pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

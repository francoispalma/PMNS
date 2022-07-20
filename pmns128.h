#ifndef MPPMNS_H
#define MPPMNS_H

#include <stdint.h>
#include <stdlib.h>

#include "structs.h"
#include "params128.h"

extern void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B);
extern void m1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A);
extern void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
extern void amns128_montg_mult_pre(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
extern void amns128_montg_mult_hyb(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);

#endif

#ifndef MONTGOM_H
#define MONTGOM_H

#include <stdint.h>
#include <stdlib.h>

#include "structs.h"
#include "params.h"

extern void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
extern void amns_montg_mult_pre(restrict poly res, const restrict poly A,
	const restrict poly B);
extern void mns_montg_int_red(restrict poly res, __int128* R);

#endif

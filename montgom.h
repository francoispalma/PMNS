#ifndef MONTGOM_H
#define MONTGOM_H

#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>

typedef struct
{
	uint16_t deg;
	int64_t* t;
} _poly, *poly;

#include "params.h"

extern void init_poly(const uint16_t deg, restrict poly* P);
void init_polys(const uint16_t deg, restrict poly* P, ...);
extern void free_poly(restrict poly P);
void free_polys(restrict poly P, ...);
void set_val(restrict poly P, int64_t val, ...);
extern void print(const restrict poly P);
extern void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);
extern void mns_montg_int_red(restrict poly res, __int128* R);

#endif

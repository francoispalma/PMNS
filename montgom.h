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

extern void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif

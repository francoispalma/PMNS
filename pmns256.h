#ifndef PMNS256_H
#define PMNS256_H

#include <stdint.h>
#include <stdlib.h>

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

#include "params256.h"

typedef struct _poly
{
	uint64_t lo[N];
	uint64_t midlo[N];
	uint64_t midhi[N];
	int64_t hi[N];
} stpoly256, *poly256;

#endif

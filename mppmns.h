#ifndef MPPMNS_H
#define MPPMNS_H

#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>

typedef struct
{
	uint16_t deg;
	int64_t* hi;
	uint64_t* lo;
} _poly128, *poly128;

#endif

/*
0|A1|0|A2
0|B1|0|B2
A1B1h|A1B1l + A1B2h + A2B1h|A1B2l + A2B1l + A2B2h|A2B2l

A1B2l + A2B1l + A2B2h|A2B2l
Q1|Q2
*/

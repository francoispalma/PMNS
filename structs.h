#ifndef AMNS_PMNS_STRUCTS_H_INCLUDED
#define AMNS_PMNS_STRUCTS_H_INCLUDED

typedef struct
{
	uint16_t deg;
	int64_t* t;
} _poly, *poly;

typedef struct
{
	uint16_t deg;
	int64_t* hi;
	uint64_t* lo;
} _poly128, *poly128;

#endif

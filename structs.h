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

extern void init_poly(const uint16_t deg, restrict poly* P);
void init_polys(const uint16_t deg, restrict poly* P, ...);
extern void free_poly(restrict poly P);
void free_polys(restrict poly P, ...);
extern void init_poly128(const uint16_t deg, restrict poly128* P);
void set_val(restrict poly P, int64_t val, ...);
extern void print(const restrict poly P);
void init_poly128s(const uint16_t deg, restrict poly128* P, ...);
extern void free_poly128(restrict poly128 P);
void free_poly128s(restrict poly128 P, ...);
void p128_set_val(restrict poly128 P, __int128 val, ...);
extern void p128_print(const restrict poly128 P);

#endif

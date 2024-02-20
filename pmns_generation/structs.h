#ifndef AMNS_PMNS_STRUCTS_H_INCLUDED
#define AMNS_PMNS_STRUCTS_H_INCLUDED

typedef struct
{
	uint16_t deg;
	int64_t* t;
} _poly, *poly;

typedef struct
{
	char sign;
	uint16_t deg;
	uint64_t* t;
} _mpnum, *mpnum;

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
void set_val(restrict poly P, int64_t val, ...);
extern void poly_print(const restrict poly P);
extern void poly_copy(restrict poly copy, const restrict poly original);
extern void init_mpnum(const uint16_t deg, restrict mpnum* P);
void init_mpnums(const uint16_t deg, restrict mpnum* P, ...);
extern void free_mpnum(restrict mpnum P);
void free_mpnums(restrict mpnum P, ...);
void mp_print(const restrict mpnum);
void mp_copy(restrict mpnum*, restrict const mpnum);
extern void init_poly128(const uint16_t deg, restrict poly128* P);
void init_poly128s(const uint16_t deg, restrict poly128* P, ...);
extern void free_poly128(restrict poly128 P);
void free_poly128s(restrict poly128 P, ...);
void p128_set_val(restrict poly128 P, __int128 val, ...);
extern void p128_print(const restrict poly128 P);
extern void p128_copy(restrict poly128 copy, const restrict poly128 original);
#endif


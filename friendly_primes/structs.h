#ifndef STRUCTS_H_INCLUDED
#define STRUCTS_H_INCLUDED

typedef struct
{
	uint16_t deg;
	int64_t* t;
} _poly, *poly;

extern void init_poly(const uint16_t deg, restrict poly* P);
void init_polys(const uint16_t deg, restrict poly* P, ...);
extern void free_poly(restrict poly P);
void free_polys(restrict poly P, ...);
void set_val(restrict poly P, int64_t val, ...);
extern void poly_print(const restrict poly P);
extern void poly_copy(restrict poly copy, const restrict poly original);

#endif


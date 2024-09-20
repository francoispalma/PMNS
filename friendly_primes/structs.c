#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"

inline void init_poly(const uint16_t deg, restrict poly* P)
{
	*P = malloc(sizeof(_poly));
	(*P)->deg = deg;
	(*P)->t = calloc(deg, sizeof(int64_t));
}

void init_polys(const uint16_t deg, restrict poly* P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		init_poly(deg, P);
		P = va_arg(args, poly*);
	} while(P != NULL);
	va_end(args);
}

inline void free_poly(restrict poly P)
{
	free(P->t);
	free(P);
}

void free_polys(restrict poly P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		free_poly(P);
		P = va_arg(args, poly);
	} while(P != NULL);
	va_end(args);
}

void set_val(restrict poly P, int64_t val, ...)
{
	va_list args;
	
	va_start(args, val);
	for(int16_t i = 0; i < P->deg; i++)
	{
		P->t[i] = val;
		val = va_arg(args, int64_t);
	}
	va_end(args);
}

inline void poly_print(const restrict poly P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[P->deg - 1]);
}

inline void poly_copy(restrict poly copy, const restrict poly original)
{
	memcpy(copy->t, original->t, original->deg * sizeof(int64_t));
}

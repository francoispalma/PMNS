#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>

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

inline void print(const restrict poly P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[P->deg - 1]);
}

inline void init_poly128(const uint16_t deg, restrict poly128* P)
{
	*P = malloc(sizeof(_poly128));
	(*P)->deg = deg;
	(*P)->hi = calloc(deg, sizeof(int64_t));
	(*P)->lo = calloc(deg, sizeof(int64_t));
}

void init_poly128s(const uint16_t deg, restrict poly128* P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		init_poly128(deg, P);
		P = va_arg(args, poly128*);
	} while(P != NULL);
	va_end(args);
}

inline void free_poly128(restrict poly128 P)
{
	free(P->hi);
	free(P->lo);
	free(P);
}

void free_poly128s(restrict poly128 P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		free_poly128(P);
		P = va_arg(args, poly128);
	} while(P != NULL);
	va_end(args);
}

void p128_set_val(restrict poly128 P, __int128 val, ...)
{
	va_list args;
	
	va_start(args, val);
	for(int16_t i = 0; i < P->deg; i++)
	{
		P->lo[i] = val;
		P->hi[i] = val >> 64;
		val = va_arg(args, __int128);
	}
	va_end(args);
}

inline void p128_print(const restrict poly128 P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%lx%lx, ", P->hi[i], P->lo[i]);
	printf("%lx%lx]\n", P->hi[P->deg - 1], P->lo[P->deg - 1]);
}

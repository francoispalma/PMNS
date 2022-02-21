#include <stdio.h>

#include "montgom.h"

static inline void init_poly(const uint16_t deg, restrict poly* P)
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

static inline void free_poly(restrict poly P)
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

static inline void print(const restrict poly P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[P->deg - 1]);
}

static inline void mns_mod_mult_ext_red(__int128* R, const restrict poly A,
	const restrict poly B, const int64_t lambda)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register const uint16_t n = A->deg;
	register uint16_t i, j;
	
	for(i = 0; i < n; i++)
	{
		for(j = 1; j < n - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[n - j];
		
		R[i] = R[i] * lambda;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B, const restrict poly M, const restrict poly M1, 
	const register int64_t lambda)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res;
	
	register uint16_t n = A->deg;
	
	__int128* R = calloc(n, sizeof(__int128));
	uint64_t* V = malloc(n * sizeof(int64_t));
	uint64_t* V2 = malloc(n * sizeof(int64_t));
	uint64_t* T = malloc(n * sizeof(int64_t));
	uint64_t* T2 = malloc(n *  sizeof(int64_t));
	
	mns_mod_mult_ext_red(R, A, B, lambda);
	
	for(int i = 0; i < n; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	mns_mod_mult_ext_red(R, res, M1, lambda);
	
	for(int i = 0; i < n; i++)
	{
		res->t[i] = R[i];
		R[i] = 0;
	}
	
	mns_mod_mult_ext_red(R, res, M, lambda);
	
	for(int i = 0; i < n; i++)
	{
		T[i] = R[i];
		T2[i] = (R[i] >> 64);
		
		T[i] = V[i] + T[i];
		res->t[i] = V2[i] + T2[i] + (T[i] < V[i]); 
	}
	
	free(R);
	free(V);
	free(V2);
	free(T);
	free(T2);
}

void __init_tests__()
{
	const register unsigned short n = 5;
	const register unsigned long lambda = 2;
	poly a, b, c, c_check, M, M1;
	init_polys(n, &a, &b, &c, &c_check, &M, &M1, NULL);
	
	set_val(a, 3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128);
	set_val(b, 1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068);
	set_val(c_check, 2302327877203981, 25683149970777821, -1798382075251775, 52479742770215631, 21994577573493812);
	set_val(M, 9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940);
	set_val(M1, 7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353);
	
	amns_montg_mult(c, a, b, M, M1, lambda);
	
	print(c);
	print(c_check);
	
	free_polys(a, b, c, c_check, M, M1, NULL);
}

int main(void)
{
	__init_tests__();
	return 0;
}

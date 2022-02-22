#include <stdio.h>
#include <time.h>

#include "montgom.h"

#define RHO 61
#define N 5
#define LAMBDA 2

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
	const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
		
		R[i] = R[i] * LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B, const restrict poly M, const restrict poly M1)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	register uint16_t n = A->deg;
	
/*	__int128* R = calloc(n, sizeof(__int128));*/
/*	uint64_t* V = malloc(n * sizeof(int64_t));*/
/*	uint64_t* V2 = malloc(n * sizeof(int64_t));*/
/*	uint64_t* T = malloc(n * sizeof(int64_t));*/
/*	uint64_t* T2 = malloc(n *  sizeof(int64_t));*/
	
	__int128 R[N] = {0};
	uint64_t V[N], V2[N], T[N], T2[N];
	
	mns_mod_mult_ext_red(R, A, B);
	
	for(int i = 0; i < n; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	mns_mod_mult_ext_red(R, res, M1);
	
	for(int i = 0; i < n; i++)
	{
		res->t[i] = R[i];
		R[i] = 0;
	}
	
	mns_mod_mult_ext_red(R, res, M);
	
	for(int i = 0; i < n; i++)
	{
		T[i] = R[i];
		T2[i] = (R[i] >> 64);
		
		T[i] = V[i] + T[i];
		res->t[i] = V2[i] + T2[i] + (T[i] < V[i]); 
	}
	
/*	free(R);*/
/*	free(V);*/
/*	free(V2);*/
/*	free(T);*/
/*	free(T2);*/
}

static inline int64_t randomint64(void)
{
	return (int64_t)(((int64_t)(rand() + rand()) << 32) ^ ((int64_t)(rand() + rand())));
}

static inline int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

static inline void randpoly(poly P)
{
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64());
}

void __init_tests__(void)
{
	poly a, b, c, c_check, M, M1;
	init_polys(N, &a, &b, &c, &c_check, &M, &M1, NULL);
	
	set_val(a, 3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128);
	set_val(b, 1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068);
	set_val(c_check, 2302327877203981, 25683149970777821, -1798382075251775, 52479742770215631, 21994577573493812);
	set_val(M, 9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940);
	set_val(M1, 7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353);
	
	amns_montg_mult(c, a, b, M, M1);
	
	print(c);
	print(c_check);
	
	free_polys(a, b, c, c_check, M, M1, NULL);
}

void __proof_of_accuracy(void)
{
	time_t seed;
	poly a, b, c, M, M1;
	init_polys(N, &a, &b, &c, &M, &M1, NULL);
	srand((unsigned)time(&seed));
	
	
	set_val(M, 9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940);
	set_val(M1, 7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353);
	
	for(register int64_t i = 0; i < 100000; i++)
	{
		randpoly(a);
		randpoly(b);
		amns_montg_mult(c, a, b, M, M1);
		
		print(a);
		print(b);
		print(c);
	}
	
	free_polys(a, b, c, M, M1, NULL);
}

int main(void)
{
	//__init_tests__();
	__proof_of_accuracy();
	return 0;
}

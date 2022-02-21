#include <stdio.h>

#include "montgom.h"

void init_poly(uint16_t deg, poly* P)
{
	*P = malloc(sizeof(_poly));
	(*P)->deg = deg;
	(*P)->t = calloc(deg, sizeof(int64_t));
}

void init_polys(uint16_t deg, poly* P, ...)
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

void free_poly(poly P)
{
	free(P->t);
	free(P);
}

void free_polys(poly P, ...)
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

void set_val(poly P, int64_t val, ...)
{
	va_list args;
	
	va_start(args, val);
	for(int16_t i = P->deg - 1; i >= 0; i--)
	{
		P->t[i] = val;
		val = va_arg(args, int64_t);
	}
	va_end(args);
}

void print(poly P)
{
	printf("[");
	for(int16_t i = P->deg - 1; i > 0; i--)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[0]);
}

void mns_mod_mult_ext_red(__int128* R, poly op1, poly op2, int64_t lambda)
{
	register uint16_t n = op1->deg;
	uint16_t i, j;
	__int128 *A, *B, lam = (__int128)(lambda);
	A = malloc(n * sizeof(__int128));
	B = malloc(n * sizeof(__int128));
	
	for(i = 0; i < n; i++)
	{
		A[i] = (__int128) op1->t[i];
		B[i] = (__int128) op2->t[i];
	}
	
	for(i = 0; i < n; i++)
	{
		for(j = 1; j < n - i; j++)
		{
			R[i] += A[i + j] * B[n - j];
			printf("%d %d 0x%lx%016lx\n", i, j, (R[i] >> 64), R[i]);
		}
		
		printf("Before lambda 0x%lx%016lx\n", (R[i] >> 64), ((R[i]<<64)>>64));
		R[i] = (__int128)((__int128)(R[i]) * (__int128)(lam));
		printf("After lambda 0x%lx%016lx\n", (R[i] >> 64), R[i]);
		
		for(j = 0; j < i + 1; j++)
			R[i] += A[j] * B[i - j];
	}
	
	free(A);
	free(B);
}

void amns_montg_mult(poly res, poly A, poly B, poly M, poly M1, int64_t lambda)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res;
	
	register uint16_t n = A->deg;
	
	__int128* R = calloc(n, sizeof(__int128));
	uint64_t* Q = calloc(n, sizeof(int64_t));
	uint64_t* Q2 = calloc(n, sizeof(int64_t));
	
	mns_mod_mult_ext_red(R, A, B, lambda);
	
	for(int i = 0; i < n; i++)
	{
		Q[i] = R[i];
		Q2[i] = (R[i] >> 64);
		printf("0x%lx%016lx\n", Q2[i], Q[i]);
	}
	
	free(R);
	free(Q);
	free(Q2);
}

void __init_tests__()
{
	unsigned short n = 5;
	unsigned long lambda = 2;
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

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

void amns_montg_mult(poly res, poly op1, poly op2, poly M, poly M1)
{
	// Function that multiplies op1 by op2 using the montgomery approach in an
	// amns. Puts the result in res;
	
	register uint16_t n = op1->deg;
	__int128* R = calloc(n, sizeof(__int128));
	
	free(R);
}

void __init_tests__()
{
	unsigned short n = 5;
	poly a, b, c, c_check, M, M1;
	init_polys(n, &a, &b, &c, &c_check, &M, &M1, NULL);
	
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

int main(void)
{
	__init_tests__();
	return 0;
}

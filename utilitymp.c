// This module is to provide multiprecision operations without relying on gmp
// This allows the code to compile even on gmp-less environments.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "structs.h"

void __print128(register const __int128 Val)
{
	int64_t hi = Val >> 64;
	uint64_t lo = Val;
	printf("0x%lx%016lx\n", hi, lo);
}

void mp_print(const restrict poly P)
{
	printf("%lx", P->t[P->deg - 1]);
	for(register uint16_t i = 1; i < P->deg; i++)
		printf("%016lx", P->t[P->deg - 1 - i]);
	printf("\n"); 
}

void mp_reduce(restrict poly A)
{
	while(A->deg > 1 && A->t[A->deg - 1] == 0) --A->deg;
}

void convert_string_to_multipre(restrict poly* res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// multiprecision number. Despite the use of the poly structure, it is not
	// a polynomial and is not in a PMNS system.
	
	uint16_t tabsize = 0, j;
	int16_t i = 0;
	char store[17];
	
	store[16] = '\0';
	
	if(string[0] != '\0' && string[1] != '\0' && string[0] == '0' && string[1] == 'x')
		string += 2;
	
	while(string[0] == '0')
		string++;
	
	while(string[i] != '\0')
	{
		if(i % 16 == 0) ++tabsize;
		++i;
	}
	if (tabsize > (*res)->deg)
	{
		free_poly(*res);
		init_poly(tabsize, res);
	}
	(*res)->deg = tabsize;
	
	// we keep i from the end of the last loop on purpose.
	for(i = i - 1; i > -1; i -= 16)
	{
		for(j = 0; j < 16; j++)
		{
			if(i - 15 + j >= 0)
				store[j] = string[i - 15 + j];
			else
				store[j] = '0';
		}
		(*res)->t[tabsize - 1 - (i / 16)] = strtoul(store, NULL, 16);
	}
	
	mp_reduce(*res);
}

void mp_copy(restrict poly* A, restrict const poly B)
{
	// Copy B into A.
	if((*A)->deg < B->deg)
	{
		free_poly(*A);
		init_poly(B->deg, A);
	}
	
	// In case of deg(A) > deg(B)
	(*A)->deg = B->deg;
	
	// Copy values in each limb.
	for(register uint16_t i = 0; i < B->deg; i++)
		(*A)->t[i] = B->t[i];
}

int8_t mp_ucomp(restrict poly A, restrict poly B)
{
	// compares A and B. Returns 0 if equal, 1 if A > B and -1 if B < A.
	// doesn't check sign.
	
	register int32_t i;
	
	mp_reduce(A);
	mp_reduce(B);

	if(A->deg > B->deg)
		return 1;
	else
	{
		if(B->deg > A->deg)
			return -1;
		else
		{			
			for(i = A->deg - 1; i > -1; i--)
			{
				if(((uint64_t)A->t[i]) > ((uint64_t)B->t[i]))
					return 1;
				if(((uint64_t)A->t[i]) < ((uint64_t)B->t[i]))
					return -1;
			}
			return 0;
		}
	}
}

int8_t mp_comp(restrict poly A, restrict poly B)
{
	// compares A and B. Returns 0 if equal, 1 if A > B and -1 if B < A.
	// checks sign.
	
	mp_reduce(A);
	mp_reduce(B);

	if((A->t[A->deg - 1] & 0x8000000000000000) &&
			(!(B->t[B->deg - 1] & 0x8000000000000000)))
		return -1;
	if((!(A->t[A->deg - 1] & 0x8000000000000000)) &&
			(B->t[B->deg - 1] & 0x8000000000000000))
		return 1;

	return mp_ucomp(A, B);
}

void mp_leftshift(restrict poly* A)
{
	// Shifts *A to the left once.
	
	poly aux;
	init_poly(0, &aux);
	
	mp_reduce(*A);
	mp_copy(&aux, *A);
	
	if((*A)->t[(*A)->deg - 1] & 0x8000000000000000)
	{
		free_poly(*A);
		init_poly(aux->deg + 1, A);
		(*A)->t[(*A)->deg - 1] = 1;
	}
	
	for(register uint16_t i = 0; i < aux->deg - 1; i++)
		(*A)->t[aux->deg - 1 - i] = (aux->t[aux->deg - 1 - i] << 1) +
			((aux->t[aux->deg - 2 - i] & 0x8000000000000000) != 0);
	
	(*A)->t[0] = aux->t[0] << 1;
	
	free_poly(aux);
}

void mp_rightshift(restrict poly A)
{
	// Shifts A to the right once.
	
	mp_reduce(A);
	
	for(register uint16_t i = 0; i < A->deg - 1; i++)
		A->t[i] = (((uint64_t)A->t[i]) >> 1) | ((A->t[i + 1] & 1) << 63);
	A->t[A->deg - 1] = (((uint64_t)A->t[A->deg - 1]) >> 1);
	
	mp_reduce(A);
}

void mp_alignleft(restrict poly* A, uint16_t deg)
{
	poly aux;
	
	mp_reduce(*A);
	init_poly((*A)->deg, &aux);
	
	if((*A)->deg < deg)
	{
		mp_copy(&aux, *A);
		free_poly(*A);
		init_poly(deg, A);
		for(register uint16_t i = 0; i < aux->deg; i++)
			(*A)->t[i + deg - aux->deg] = aux->t[i];
	}
	
	free_poly(aux);
}

void mp_add(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 + op2 in res.

	const uint16_t MAXDEG = (op1->deg < op2->deg ? op2->deg : op1->deg) + 1; // -
		//((op1->t[op1->deg - 1] < 0) || (op2->t[op2->deg - 1] < 0));
	register uint16_t i, j;
	uint64_t stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	(*res)->deg = MAXDEG;
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = ((uint64_t) (*res)->t[i]);
		(*res)->t[i] += op2->t[i];
		
		j = i;
		while(stok > ((uint64_t) (*res)->t[j]) && j < MAXDEG - 1)
		{
			++j;
			stok = ((uint64_t) (*res)->t[j]);
			(*res)->t[j] = ((uint64_t) (*res)->t[j]) + 1;
		}
	}
	
	mp_reduce(*res);
}

void mp_sub(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 - op2 in res.

	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	register uint16_t i, j;
	uint64_t stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	(*res)->deg = MAXDEG;
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = ((uint64_t) (*res)->t[i]);
		(*res)->t[i] -= op2->t[i];
		
		j = i;
		while(stok < ((uint64_t) (*res)->t[j]) && j < MAXDEG - 1)
		{
			++j;
			stok = ((uint64_t) (*res)->t[j]);
			(*res)->t[j] = ((uint64_t) (*res)->t[j]) - 1;
		}
	}
	
	mp_reduce(*res);
}

void mp_usub(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 - |op2| in res.

	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	register uint16_t i, j;
	const int8_t sign = op2->t[op2->deg - 1] >= 0 ? 1 : -1;
	uint64_t stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	(*res)->deg = MAXDEG;
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = ((uint64_t) (*res)->t[i]);
		(*res)->t[i] -= sign * op2->t[i];
		
		j = i;
		while(stok < ((uint64_t) (*res)->t[j]) && j < MAXDEG - 1)
		{
			++j;
			stok = ((uint64_t) (*res)->t[j]);
			(*res)->t[j] = ((uint64_t) (*res)->t[j]) - 1;
		}
	}
	
	mp_reduce(*res);
}

void mp_mult(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 * op2 into res.

	register uint16_t i, j, k;
	unsigned __int128 R[op1->deg + op2->deg], stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < op1->deg + op2->deg)
	{
		free_poly(*res);
		init_poly(op1->deg + op2->deg, res);
	}
	(*res)->deg = op1->deg + op2->deg;
	
	for(i = 0; i < op1->deg + op2->deg; i++)
		R[i] = 0;

	for(i = 0; i < op1->deg; i++)
	{
		for(j = 0; j < op2->deg; j++)
		{
			k = i + j;
			stok = R[k];
			R[k] += (unsigned __int128) ((uint64_t) op1->t[i]) * ((uint64_t)op2->t[j]);
			while(stok > R[k])
			{
				++k;
				stok = R[k];
				R[k] += 1;
			}
		}
	}
	
	(*res)->t[0] = R[0];
	for(i = 1; i < op1->deg + op2->deg; i++)
	{
		k = i;
		stok = R[i];
		R[i] += R[i - 1] >> 64;
		while(stok > R[k])
		{
			++k;
			stok = R[k];
			R[k] += 1;
		}
		(*res)->t[i] = R[i];
	}
	
	mp_reduce(*res);
}

void mp_mod(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 % op2 into res.
	
	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	poly X, R;
	init_polys(MAXDEG, &X, &R, NULL);
	
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	
	mp_copy(&X, op2);
	mp_alignleft(&X, op1->deg);
	while(mp_ucomp(op1, X) == 1)
		mp_leftshift(&X);
	while(mp_ucomp(op1, X) == -1)
		mp_rightshift(X);
	
	if(op1->t[op1->deg - 1] > 0)
		mp_sub(res, op1, X);
	else
	{
		mp_leftshift(&X);
		mp_usub(res, X, op1);
	}
	
	while(mp_ucomp(*res, op2) == 1)
	{
		mp_copy(&R, *res);
		while(mp_ucomp(X, *res) == 1)
			mp_rightshift(X);
		mp_sub(res, R, X);
	}
	
	mp_reduce(*res);
	free_polys(X, R, NULL);
}


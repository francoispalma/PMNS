// This module is to provide multiprecision operations without relying on gmp
// This allows the code to compile even on gmp-less environments.
// The algorithms used here are not fast.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "structs.h"


void __print128(register const __int128 Val)
{
	int64_t hi = Val >> 64;
	uint64_t lo = Val;
	printf("%lx%016lx", hi, lo);
}

static inline void mp_reduce(mpnum A)
{
	while(A->deg > 1 && A->t[A->deg - 1] == 0) --A->deg;
}

void convert_string_to_binary(mpnum* res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// multiprecision number.
	
	uint16_t tabsize = 0, j;
	int16_t i = 0;
	char store[17], sign = 1;
	
	store[16] = '\0';
	
	// We deal with potential negative numbers.
	if(string[0] == '-')
	{
		sign = -1;
		string++;
	}
	
	// We accept 0x in front but will treat it as hex anyway so we remove it.
	if(string[0] != '\0' && string[1] != '\0' && string[0] == '0' && string[1] == 'x')
		string += 2;
	
	// We remove any 0 in front of it.
	while(string[0] == '0')
		string++;
	
	// We check how many limbs will be needed.
	while(string[i] != '\0')
	{
		if(i % 16 == 0) ++tabsize;
		++i;
	}
	
	// We reallocate if needed.
	if (tabsize > (*res)->deg)
	{
		free_mpnum(*res);
		init_mpnum(tabsize, res);
	}
	(*res)->deg = tabsize;
	(*res)->sign = sign;
	
	// We keep i from the end of the last loop on purpose for the fill order.
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

_Bool mp_iszero(const mpnum A)
{
	// Returns 1 if A = 0 else returns 0.
	mp_reduce(A);
	
	if(A->deg == 1 && A->t[0] == 0)
		return 1;
	else
		return 0;
}

int8_t mp_ucomp(mpnum A, mpnum B)
{
	// Compares |A| and |B|. Returns 0 if equal, 1 if |A| > |B| and -1 if 
	// |B| < |A|.
	
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
				if(A->t[i] > B->t[i])
					return 1;
				if(A->t[i] < B->t[i])
					return -1;
			}
			return 0;
		}
	}
}

int8_t mp_comp(mpnum A, mpnum B)
{
	// compares A and B. Returns 0 if equal, 1 if A > B and -1 if B < A.
	
	mp_reduce(A);
	mp_reduce(B);
	
	if(A->sign == B->sign)
		return mp_ucomp(A, B) * A->sign;
	else
		return A->sign;
}

void mp_leftshift(mpnum* A)
{
	// Shifts *A to the left once.
	
	mpnum aux;
	init_mpnum(0, &aux);
	
	mp_reduce(*A);
	mp_copy(&aux, *A);
	
	if((*A)->t[(*A)->deg - 1] & 0x8000000000000000)
	{
		free_mpnum(*A);
		init_mpnum(aux->deg + 1, A);
		(*A)->t[(*A)->deg - 1] = 1;
		(*A)->sign = aux->sign;
	}
	
	for(register uint16_t i = 0; i < aux->deg - 1; i++)
		(*A)->t[aux->deg - 1 - i] = (aux->t[aux->deg - 1 - i] << 1) +
			((aux->t[aux->deg - 2 - i] & 0x8000000000000000) != 0);
	
	(*A)->t[0] = aux->t[0] << 1;
	
	free_mpnum(aux);
}

void mp_rightshift(mpnum A)
{
	// Shifts A to the right once.
	
	mp_reduce(A);
	
	for(register uint16_t i = 0; i < A->deg - 1; i++)
		A->t[i] = (((uint64_t)A->t[i]) >> 1) | ((A->t[i + 1] & 1) << 63);
	A->t[A->deg - 1] = (((uint64_t)A->t[A->deg - 1]) >> 1);
	
	mp_reduce(A);
}

void mp_alignleft(mpnum* A, uint16_t deg)
{
	mpnum aux;
	
	mp_reduce(*A);
	init_mpnum((*A)->deg, &aux);
	
	if((*A)->deg < deg)
	{
		mp_copy(&aux, *A);
		free_mpnum(*A);
		init_mpnum(deg, A);
		(*A)->sign = aux->sign;
		for(register uint16_t i = 0; i < aux->deg; i++)
			(*A)->t[i + deg - aux->deg] = aux->t[i];
	}
	
	free_mpnum(aux);
}

void mp_uadd(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of |op1| + |op2| in res.
	
	const uint16_t MAXDEG = (op1->deg < op2->deg ? op2->deg : op1->deg) + 1;
	register uint16_t i, j;
	uint64_t stok;
	
	// We put res at 0 with the right size.
	free_mpnum(*res);
	init_mpnum(MAXDEG, res);
	(*res)->sign = 1;
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = (*res)->t[i];
		(*res)->t[i] += op2->t[i];
		
		// We deal with the carry.
		j = i;
		while(j < MAXDEG - 1 && stok > (*res)->t[j])
		{
			++j;
			stok = (*res)->t[j];
			(*res)->t[j] += 1;
		}
	}
	
	mp_reduce(*res);
}

void mp_usub(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of |op1| - |op2| in res.
	
	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	register uint16_t i, j;
	uint64_t stok;
	
	// We check if the degree is high enough. If it isn't we fix the problem.
	free_mpnum(*res);
	init_mpnum(MAXDEG, res);
	(*res)->sign = mp_ucomp(op1, op2);
	(*res)->sign += ((*res)->sign == 0);
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = (*res)->t[i];
		(*res)->t[i] -= op2->t[i];
		
		j = i;
		while(j < MAXDEG - 1 && stok < (*res)->t[j])
		{
			++j;
			stok = (*res)->t[j];
			(*res)->t[j] -= 1;
		}
	}
	
	mp_reduce(*res);
}

void mp_add(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of op1 + op2 in res.
	
	if(op1->sign == op2->sign)
	{
		mp_uadd(res, op1, op2);
		(*res)->sign = op1->sign;
	}
	else
	{
		if(mp_ucomp(op1, op2) == 1)
		{
			mp_usub(res, op1, op2);
			(*res)->sign *= op1->sign;
		}
		else
		{
			mp_usub(res, op2, op1);
			(*res)->sign *= op2->sign;
		}
	}
}

void mp_sub(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of op1 - op2 in res.
	
	op2->sign = -op2->sign;
	mp_add(res, op1, op2);
	op2->sign = -op2->sign;
}

void mp_mult(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of op1 * op2 into res.
	
	register uint16_t i, j, k;
	unsigned __int128 R[op1->deg + op2->deg], stok;
	_Bool carry;
	
	// We put res at 0 with the right size.
	free_mpnum(*res);
	init_mpnum(op1->deg + op2->deg, res);
	(*res)->sign = op1->sign * op2->sign;
	
	for(i = 0; i < op1->deg + op2->deg; i++)
		R[i] = 0;

	for(i = 0; i < op1->deg; i++)
	{
		for(j = 0; j < op2->deg; j++)
		{
			k = i + j;
			carry = __builtin_add_overflow((unsigned __int128) op1->t[i] * op2->t[j],
					R[k], &stok);
			while(carry && k + 2 < (*res)->deg )
			{
				R[k] = stok;
				k++;++k; // We add the carry at k+2 because of the reconstruction later.
				carry = __builtin_add_overflow(R[k], (unsigned __int128) 1, &stok);
			}
			R[k] = stok;
		}
	}
	
	(*res)->t[0] = R[0];
	for(i = 1; i < op1->deg + op2->deg; i++)
	{
		k = i;
		carry = __builtin_add_overflow((unsigned __int128) (R[k - 1] >> 64),
					R[k], &stok);
		while(carry && k < (*res)->deg)
		{
			R[k] = stok;
			++k;
			carry = __builtin_add_overflow((unsigned __int128) 1, R[k], &stok);
		}
		R[k] = stok;
		(*res)->t[i] = R[i];
	}
	
	mp_reduce(*res);
}

void mp_mod(mpnum* res, const mpnum op1, const mpnum op2)
{
	// Puts the result of op1 % op2 into res. Note that this is not the fastest
	// modular algorithm possible for now but this is not the focus.
	
	if(op2->sign == -1)  // modulus of a negative number makes no sense
		return;
	
	mpnum X, R;
	
	// We put res at 0 with the right size (res is always positive).
	free_mpnum(*res);
	init_mpnum(op2->deg, res);
	
	if(mp_ucomp(op1, op2) == -1)
	{
		if(op1->sign == -1)
		{
			mp_add(res, op1, op2);
			return;
		}
		else
		{
			mp_copy(res, op1);
			return;
		}
	}
	
	init_mpnums(op1->deg, &X, &R, NULL);
	
	// We proceed using the Russian Peasant algorithm.
	
	// We put X such that the highest bit of X matches the highest bit of op1.
	// To do that we first align the degrees.
	mp_copy(&X, op2);
	mp_alignleft(&X, op1->deg);
	// Then we leave it such that X < |op1| < 2X
	while(mp_ucomp(op1, X) == 1)
		mp_leftshift(&X);
	while(mp_ucomp(op1, X) == -1)
		mp_rightshift(X);
	if(__builtin_clzll(X->t[op1->deg - 1]) != __builtin_clzll(op1->t[op1->deg - 1]))
		mp_leftshift(&X);
	
	// If op1 is negative we want to subtract -X
	X->sign = op1->sign;
	mp_sub(res, op1, X);
	
	// Then we continue.
	while(mp_ucomp(*res, op2) == 1)
	{
		mp_copy(&R, *res);
		while(mp_ucomp(X, *res) == 1)
			mp_rightshift(X);
		if(__builtin_clzll(X->t[X->deg - 1]) != __builtin_clzll((*res)->t[X->deg - 1]))
			mp_leftshift(&X);
		X->sign = (*res)->sign;
		mp_sub(res, R, X);
	}
	
	// At the very end if we have a negative value we add op2 one last time.
	if((*res)->sign == -1)
	{
		mp_copy(&R, *res);
		mp_add(res, R, op2);
	}
	
	mp_reduce(*res);
	free_mpnums(X, R, NULL);
}

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

void random_mp(mpnum ret)
{
	for(int i = 0; i < ret->deg; i++)
		ret->t[i] = randomint64();
}


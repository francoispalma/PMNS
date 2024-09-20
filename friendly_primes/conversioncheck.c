#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include <time.h>

#include "hpmns.h"
#include "params.h"

#define NSAMPLES 1000

int64_t randomint64(void)
{
	// Function that generates a random integer on 64 bits.
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

void randhlimbs(uint64_t g[], const uint16_t limbs)
{
	// Function that generates a random limb array.
	for(uint64_t i = 0; i < limbs; i++)
		g[i] = randomint64();
}

static inline int8_t _sign(int64_t var)
{
	return (var >= 0) - (var < 0);
}

void horner_modulo(mp_limb_t res[], const poly P, const uint16_t plimbs)
{
	// Function that evaluates the polynomial P in gamma.
	// Uses the horner algorithm.
	// Reduces the result modulo __P__ using Montgomery CIOS.
	
	int8_t sign = _sign(P->t[N-1]);
	mp_limb_t A[N+1] = {sign*P->t[N-1], 0};
	uint16_t sizeA = 1;
	mp_limb_t carry;
	
	// Horner evaluation.
	for(int i = N - 2; i >= 0; i--)
	{
		carry = mpn_mul_1(A, A, sizeA, GAMMA);
		A[sizeA] = carry;
		sizeA += (carry != 0);
		if(_sign(P->t[i]) == sign)
			mpn_add_1(A, A, sizeA + 1, sign*P->t[i]);
		else
		{
			carry = mpn_sub_1(A, A, sizeA + 1, -sign*P->t[i]);
			if(carry)
			{
				sign = -sign;
				mpn_mul_1(A, A, sizeA + 1, -1);
				sizeA = 1;
			}
		}
	}
	
	// We multiply by alpha.
	A[sizeA] = mpn_mul_1(A, A, sizeA, ALPHA);
	
	// Montgomery reduction.
	mp_limb_t Q[N + 1] = {0};
	Q[N] = mpn_mul_1(Q, __P__, N, __PINVERSE__ * A[0]);
	mpn_add_n(Q, Q, A, N + 1);
	
	if(sign == 1)
		for(int i = 0; i < plimbs; i++)
			res[i] = Q[i + 1];
	else
		mpn_sub_n(res, __P__, Q+1, plimbs);
}

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	int64_t tabs[NSAMPLES][N] = {{0}};
	_poly dummy[NSAMPLES];
	for(int i = 0; i < NSAMPLES; i++)
	{
		dummy[i].t = tabs[i];
	}
	
	uint16_t plimbs = N;
	while (__P__[plimbs - 1] == 0)
		plimbs -= 1;
	
	uint64_t toconv[NSAMPLES][N] = {{0}};
	uint64_t convback[NSAMPLES][N] = {{0}};
	
	for(int i = 0; i < NSAMPLES; i++)
	{
		randhlimbs(toconv[i], plimbs);
		toconv[i][plimbs-1] = toconv[i][plimbs-1] % __P__[plimbs - 1];
	}
	
	uint64_t correct = 0;
	for(int i = 0; i < NSAMPLES; i++)
	{
		convert_binary_to_pmns(dummy + i, toconv[i]);
		horner_modulo(convback[i], dummy + i, plimbs);
		correct += (mpn_cmp(toconv[i], convback[i], plimbs) == 0);
	}
	
	printf("Correct: %ld/%d\n", correct, NSAMPLES);
	
	uint64_t morethanrho = 0;
	for(int i = 0; i < NSAMPLES; i++)
		for(int j = 0; j < N; j++)
			morethanrho += ((dummy[i].t[j] > (int64_t) RHO) || (dummy[i].t[j] < -((int64_t) RHO)));
	
	printf("More than rho: %ld\n", morethanrho);
	
	for(uint64_t i = 0; i < 10000001; i++)
	{
		printf("\b%ld\r", i);
		randhlimbs(toconv[0], plimbs);
		toconv[0][plimbs-1] = toconv[0][plimbs-1] % __P__[plimbs - 1];
		convert_binary_to_pmns(dummy, toconv[0]);
		horner_modulo(convback[0], dummy, plimbs);
		if(mpn_cmp(toconv[0], convback[0], plimbs) != 0)
		{
			printf("%ld not correct\n", i);
			break;
		}
	}
	
	printf("\n");
	
	return 0;
}

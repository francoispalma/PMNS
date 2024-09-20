#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

#include "hpmns.h"
#include "params.h"

void pmns_mod_mult_ext_red(__int128* restrict R, const restrict poly A,
		const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = alpha*X^n - lambda as a polynomial used for reduction. Result in R
	
	// Toeplitz decomposition is only worth it when N > 7 from our tests.
	#if N > 7 && (!(N % 2) || !(N % 3))
		int64_t matr[2*N - 1];
		
		for(int i = 0; i < N-1; i++)
		{
			matr[i + N - 1] = ALPHA * B->t[i];
			matr[i] = B->t[1 + i] * LAMBDA;
		}
		matr[2*N - 2] = B->t[N - 1] * ALPHA;
		toeplitz_vm(R, A->t, matr);
	#else
		// The speed gain from these is probably marginal.
		#if ALPHA != 1
			int64_t atab[N];
			for(int i = 0; i < N; i++)
				atab[i] = A->t[i]*ALPHA;
		#else
			int64_t *atab;
			atab = A->t;
		#endif
		#if LAMBDA == 1
			int64_t *btab;
			btab = B->t;
		#else
			int64_t btab[N];
			for(int i = 1; i < N; i++)
				btab[i] = B->t[i]*LAMBDA;
		#endif
		for(int i = 0; i < N; i++)
		{
			R[i] = (__int128) atab[0] * B->t[i];
			for(int j = 1; j < i + 1; j++)
				R[i] += (__int128) atab[j] * B->t[i - j];
		}
		for(int i = 0; i < N - 1; i++)
		{
			for(int j = 1; j < N - i; j++)
				R[i] += (__int128) A->t[i + j] * btab[N - j];
		}
	#endif
}


static inline void pmns_mult_by_g(__int128* restrict R, const int64_t* restrict A)
{
	// Function that multiplies A by (sparse) G. Result in R.
	
	R[0] += -(__int128)A[0] * GAMMA + (__int128)A[N - 1] * LAMED;
	for(int i = 1; i < N - 1; i++)
	{
		R[i] += A[i - 1] - (__int128)A[i] * GAMMA;
	}
	R[N - 1] += -(__int128)A[N - 1] * (GIMEL * ALEPH) + A[N - 2];
}

static inline void pmns_mult_by_g1(int64_t* restrict R, const __int128* restrict A)
{
	// Function that multiplies A by G1. Result in R.
	
	#ifdef SPARSEG1
		// If the inverse matrix is sparse we ignore the 0s.
		for(int i = 0; i < N - 2; i++)
		{
			R[i] = -(uint64_t)A[i + 1] - (uint64_t)A[i + 2] * GAMMA;
		}
		R[N - 2] = (uint64_t)A[0] * GAMMALAMM1 - (uint64_t)A[N - 1];
		R[N - 1] = (uint64_t)A[0] * ONELAMM1 + (uint64_t)A[1] * GAMMALAMM2;
	#else
		// If the inverse matrix isn't sparse we first compute the last coefficient.
		uint64_t Z = 0;
		for(int j = 0; j < N; j++)
			Z += (uint64_t)A[j] * (uint64_t)lastcol[j];
		R[N - 1] = Z;
		
		// Each subsequent coefficient is computed from it in a linear fashion.
		Z *= ALEPH*GIMEL;
		Z -= (uint64_t)A[N - 1];
		R[N - 2] = Z;
		for(int i = 1; i < N - 1; i++)
		{
			Z *= GAMMA;
			Z -= (uint64_t)A[N - 1 - i];
			R[N - 2 - i] = Z;
		}
	#endif
}

static inline void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Montgomery method.
	int64_t T[N];
	
	// T <- R times G
	pmns_mult_by_g1(T, R);
	
	// R <- R + T times G1
	pmns_mult_by_g(R, T);
	
	// res <- R divided by phi.
	for(int i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64);
}

void pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// pmns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	// R <- A times B mod E
	pmns_mod_mult_ext_red(R, A, B);
	
	// res <- Gmont-like(R)
	pmns_montg_int_red(res, R);
}

void convert_binary_to_pmns(restrict poly res, const uint64_t* restrict op)
{
	// Function that converts an integer into a polynomial in our representation
	// system. Assumes phi = 2^64.
	uint8_t counter;
	uint16_t curr;
	const uint64_t theta = (1ULL<<THETALOG2);
	__int128 R[N] = {0};
	
	// Theta-radix decomposition of op
	curr = 0;
	counter = 0;
	for(uint16_t i = 0; i < N; i++)
	{
		// If we don't add '* (counter != 0)' the result is incorrect sometimes.
		res->t[i] = ((op[curr]>>counter) | ((op[curr + 1] << (64-counter)) * (counter != 0))) & (theta - 1);
		counter += THETALOG2;
		if(counter >= 64)
		{
			counter -= 64;
			curr += 1;
		}
	}
	
	// Multiplication by precomputed polynomials representing powers of Theta.
	for(uint16_t i = 0; i < N; i++)
		for(uint16_t j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	// res <- Gmont-like(R)
	pmns_montg_int_red(res, R);
}

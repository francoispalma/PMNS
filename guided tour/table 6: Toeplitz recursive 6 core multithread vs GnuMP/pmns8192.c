#include <stdint.h>
#include <omp.h>

#include "pmns8192.h"
#include "params8192.h"

void toeplitz_vm_63x63(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 63
	for(int i = 0; i < 63; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 63
		for(int j = 0; j < 63; j++)
			rop[i] += (__int128) vect[j] * matr[63 - 1 - j + i];
	}
}

void toeplitz_vm3x3_189x189(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[63], t1[63], t2[63], t3[63], t4[63], t5[63];
	
	int64_t m03, v1m2[63], v0m2[63], v0m1[63], m034[125], m013[125], m012[125], m0[125], m1[125], m3[125];
	
	
	#pragma GCC unroll 63
	for(int i = 0; i < 63; i++)
	{
		v1m2[i] = vect[63 + i] - vect[126 + i];
		v0m1[i] = vect[i] - vect[63 + i];
		v0m2[i] = vect[i] - vect[126 + i];
	}
	
	#pragma GCC unroll 125
	for(int i = 0; i < 125; i++)
	{
		m0[i] = matr[i + 126];
		m1[i] = matr[i + 189];
		m3[i] = matr[i + 63];
		m03 = m0[i] + m3[i];
		m034[i] = m03 + matr[i];
		m013[i] = m03 + m1[i];
		m012[i] = m0[i] + m1[i] + matr[i + 252];
	}
	
	#pragma omp parallel sections num_threads(6)
	{
		#pragma omp section
		toeplitz_vm_63x63(t0, vect + 126, m034);
		#pragma omp section
		toeplitz_vm_63x63(t1, vect + 63, m013);
		#pragma omp section
		toeplitz_vm_63x63(t2, vect, m012);
		#pragma omp section
		toeplitz_vm_63x63(t3, v1m2, m3);
		#pragma omp section
		toeplitz_vm_63x63(t4, v0m2, m0);
		#pragma omp section
		toeplitz_vm_63x63(t5, v0m1, m1);
	}
	
	#pragma GCC unroll 192
	for(int i = 0; i < 63; i++)
	{
		rop[i] = t0[i] + t3[i] + t4[i];
		rop[i + 63] = t1[i] - t3[i] + t5[i];
		rop[i + 126] = t2[i] - t4[i] - t5[i];
	}
}

static inline void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	int64_t matrix[2*N - 1];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];
	
	toeplitz_vm3x3_189x189(R, A->t, matrix);
}

static inline void m_pmns_mult(__int128* restrict R, int64_t* restrict A)
{
	// Multiplies A and M's companion matrix, result in R.
	
	__int128 somme;
	
	
	#pragma omp parallel for shared(A, matrM, R) private(somme) num_threads(6)
	for(int i = 0; i < N; i++)
	{
		somme = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			somme += (__int128) A[j] * matrM[N - 1 - j + i];
		
		R[i] += somme;
	}
}

static inline void m1_pmns_mult(int64_t* restrict R, __int128* restrict A)
{
	// Multiplies A and M^-1's companion matrix mod PHI, result in R.
	
	#pragma omp parallel for shared(A, matrM1, R) num_threads(6)
	for(int i = 0; i < N; i++)
	{
		R[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			R[i] += (int64_t) A[j] * matrM1[N - 1 - j + i];
	}
}

static inline void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Polynomial Montgomery-like method.
	int64_t T[N];
	register uint16_t i;
	
	m1_pmns_mult(T, R);
	
	m_pmns_mult(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

static inline void pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the polynomial reduction in a
	// pmns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	pmns_mod_mult_ext_red(R, A, B);
	
	pmns_montg_int_red(res, R);
}

void pmns8192_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// 8192-bit modular multiplication using the polynomial reduction algorithm.
	pmns_montg_mult(res, A, B);
}


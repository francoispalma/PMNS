#include <stdint.h>
#include <omp.h>

static inline void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	int64_t matrix[2*N - 1];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		somme = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			somme += (__int128) A->t[j] * matrix[N - 1 - j + i];
		R[i] = somme;
	}
}

static inline void b_pmns_mult(__int128* restrict R,
	int64_t* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	
	__int128 somme;
	
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		somme = (__int128)B[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (__int128)B[j][i]*A[j];
		
		R[i] += somme;
	}
}

static inline void b1_pmns_mult(int64_t* restrict R,
	__int128* restrict A)
{
	// Vector-Matrix multiplication between A and B^-1, result in R.
	
	uint64_t somme;
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		somme = (uint64_t)B1[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (uint64_t)B1[j][i]*A[j];
		R[i] = somme;
	}
}

static inline void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Polynomial Montgomery-like method.
	int64_t T[N];
	register uint16_t i;
	
	b1_pmns_mult(T, R);
	
	b_pmns_mult(R, T);
	
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

static inline void ppmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B, const uint8_t NBTHREADZ)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	int64_t matrix[2*N - 1];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];
	
	#pragma omp parallel for shared(A, matrix, R) private(somme) num_threads(NBTHREADZ)
	for(int i = 0; i < N; i++)
	{
		somme = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			somme += (__int128) A->t[j] * matrix[N - 1 - j + i];
		R[i] = somme;
	}
}

static inline void pb_pmns_mult(__int128* restrict R,
	int64_t* restrict A, const uint8_t NBTHREADZ)
{
	// Vector-Matrix multiplication between A and B, result in R.
	
	__int128 somme;
	
	
	#pragma omp parallel for shared(A, B, R) private(somme) num_threads(NBTHREADZ)
	for(int i = 0; i < N; i++)
	{
		somme = (__int128)B[0][i]*A[0];
		_PRAGMAGCCUNROLLLOOP_
		for (int j = 1; j < N; j++)
			somme += (__int128)B[j][i]*A[j];
		
		R[i] += somme;
	}
}

static inline void pb1_pmns_mult(int64_t* restrict R,
	__int128* restrict A, const uint8_t NBTHREADZ)
{
	// Vector-Matrix multiplication between A and B^-1, result in R.
	
	uint64_t somme;
	
	#pragma omp parallel for shared(A, B1, R) private(somme) num_threads(NBTHREADZ)
	for(int i = 0; i < N; i++)
	{
		somme = (uint64_t)B1[0][i]*A[0];
		_PRAGMAGCCUNROLLLOOP_
		for (int j = 1; j < N; j++)
			somme += (uint64_t)B1[j][i]*A[j];
		R[i] = somme;
	}
}

static inline void ppmns_montg_int_red(restrict poly res, __int128* restrict R, const uint8_t NBTHREADZ)
{
	// Internal reduction of R via the Polynomial Montgomery-like method.
	int64_t T[N];
	register uint16_t i;
	
	pb1_pmns_mult(T, R, NBTHREADZ);
	
	pb_pmns_mult(R, T, NBTHREADZ);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

static inline void ppmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B, const uint8_t NBTHREADZ)
{
	// Function that multiplies A by B using the polynomial reduction in a
	// pmns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	ppmns_mod_mult_ext_red(R, A, B, NBTHREADZ);
	
	ppmns_montg_int_red(res, R, NBTHREADZ);
}

#include <stdint.h>

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

static inline _Bool add_overflow(unsigned __int128* restrict a, const unsigned __int128 b)
{
	//return __builtin_add_overflow(*a, b, a);
	const unsigned __int128 tmp = *a;
	*a += b;
	return *a < tmp;
}

static inline void pmns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	unsigned __int128 A1B0_A0B1, A0B0;
	__int128 A1B1, aux1, aux2;
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		A1B0_A0B1 = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
		{
			A1B1 = (__int128) A->hi[i + j] * B->hi[N - j];
			A0B0 = (__int128) A->lo[i + j] * B->lo[N - j];
			A1B0_A0B1 += (__int128) ((__int128) ((__int128) A->hi[i + j] * B->lo[N - j]) + ((__int128) A->lo[i + j] * B->hi[N - j]));
	
			Rhi[i] += (__int128) A1B1 + add_overflow(Rlo + i, A0B0);
		}
		
		Rhi[i] += (__int128) HIGH(A1B0_A0B1) + add_overflow(Rlo + i, (__int128) (LOW(A1B0_A0B1)) << 64);
		
		aux1 = (unsigned __int128) LOW(Rlo[i]) * LAMBDA;
		aux2 = (unsigned __int128) HI(Rlo[i]) * LAMBDA + HIGH(aux1);
		Rlo[i] = ((__int128) aux2 << 64) | LOW(aux1);
		Rhi[i] = (__int128) Rhi[i] * LAMBDA + HIGH(aux2);
		
		
		A1B0_A0B1 = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			A1B1 = (__int128) A->hi[j] * B->hi[i - j];
			A0B0 = (__int128) A->lo[j] * B->lo[i - j];
			A1B0_A0B1 += (__int128) ((__int128) ((__int128) A->hi[j] * B->lo[i - j]) + ((__int128) A->lo[j] * B->hi[i - j]));
	
			Rhi[i] += (__int128) A1B1 + add_overflow(Rlo + i, A0B0);
		}
		
		Rhi[i] += (__int128) HIGH(A1B0_A0B1) + add_overflow(Rlo + i, (__int128) (LOW(A1B0_A0B1)) << 64);
	}
}


static inline void m_pmns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	unsigned __int128 A0B1;
	__int128 A1B0;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(NBTHREADZ) private(A0B1, A1B0)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			Rhi[i] += (__int128) HIGH(A[j]) * HIGH(matrM[N - 1 - j + i]) + add_overflow(Rlo + i, ((__int128) LOW(A[j]) * LOW(matrM[N - 1 - j + i])));
		
		A0B1 = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			A0B1 += ((__int128) LOW(A[j]) * HIGH(matrM[N - 1 - j + i]));
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
		{
			A1B0 = (__int128) HIGH(A[j]) * LOW(matrM[N - 1 - j + i]);
			Rhi[i] += HIGH(A1B0) + add_overflow(Rlo + i, ((__int128)(LOW(A1B0)) << 64));
		}
		
		
		Rhi[i] += HIGH(A0B1) + add_overflow(Rlo + i, ((__int128)LOW(A0B1) << 64));
	}
}

static inline void m1_pmns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Function that multiplies A by M^-1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			Rlo[i] += (__int128) A[j] * matrM1[N - 1 - j + i];
	}
}

static inline void pmns128_montg_int_red(restrict poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho using Montgomery's internal reduction algorithm.
	unsigned __int128 V[N];
	register uint16_t i;
	
	m1_pmns128_mod_mult_ext_red(V, Rlo);
	
	m_pmns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

static inline void pmns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	pmns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	pmns128_montg_int_red(res, Rhi, Rlo);
}


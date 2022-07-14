#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <string.h>

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

#include "pmns128.h"
#include "utilitymp.h"


__inline uint64_t mulx64(uint64_t x, uint64_t y, uint64_t* hi)
{
    __asm__(
        "mulx %3, %0, %1    \n\t"
        : "=&d"(x), "=&a"(y)
        : "0"(x), "1"(y)
    );

    *hi = y;
    return x;
}

static inline void mult128(__int128* Rhi, unsigned __int128* Rlo,
	const int64_t Ahi, const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	unsigned __int128 aux, tmp;
	*Rhi = (__int128) LOW(Ahi) * LOW(Bhi);
	*Rlo = (__int128) Alo * Blo;
	aux = (__int128) Alo * LOW(Bhi);
	tmp = (__int128) *Rlo + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo > tmp);
	aux = (__int128) LOW(Ahi) * Blo;
	*Rlo = (__int128) tmp + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo < tmp);
}

static inline void umult128(unsigned __int128* Rhi, unsigned __int128* Rlo,
	const uint64_t Ahi, const uint64_t Alo, const uint64_t Bhi,
	const uint64_t Blo)
{
	unsigned __int128 aux, tmp;
	*Rhi = (__int128) LOW(Ahi) * LOW(Bhi);
	*Rlo = (__int128) Alo * Blo;
	aux = (__int128) Alo * LOW(Bhi);
	tmp = (__int128) *Rlo + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo > tmp);
	aux = (__int128) LOW(Ahi) * Blo;
	*Rlo = (__int128) tmp + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo < tmp);
}

static inline void umultadd128(unsigned __int128* Rhi, unsigned __int128* Rlo,
	const uint64_t Ahi, const uint64_t Alo, const uint64_t Bhi,
	const uint64_t Blo)
{
	unsigned __int128 auxlo, auxhi;
	
	umult128(&auxhi, &auxlo, Ahi, Alo, Bhi, Blo);
	
	*Rlo = (__int128) *Rlo + auxlo;
	*Rhi = (__int128) *Rhi + auxhi + (*Rlo < auxlo);
}

static inline void multadd128x(__int128* Rhi, unsigned __int128* Rlo,
	const int64_t Ahi, const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	unsigned __int128 auxlo;
	__int128 auxhi;
	
	mult128(&auxhi, &auxlo, Ahi, Alo, Bhi, Blo);
	
	*Rlo = (__int128) *Rlo + auxlo;
	*Rhi = (__int128) *Rhi + auxhi + (*Rlo < auxlo);
}

static inline __int128 signed_unsigned_mul128(const uint64_t A, const int64_t B)
{
	return (__int128) A * B;
}

static inline __int128 signed_unsigned_mul128x(const uint64_t A, const int64_t B)
{
	int64_t q;
	uint64_t r;
	__int128 tmp;
	tmp = (__int128) A * LOW(B);
	r = tmp;
	q = HIGH(tmp);
	q -= A;
	return (__int128) ((__int128) q << 64) | r; 
}

static inline void multadd128_old(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R.
	unsigned __int128 A0B0, tmplo;
	__int128 A1B1, A1B0, A0B1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0);
	
	tmplo = (__int128) A0B0 + ((__int128)(LOW(A0B1) + LOW(A1B0)) << 64);
	*Rhi += (__int128) aux2 + A1B1 + add_overflow(Rlo, tmplo);
}

static inline void multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R.
	unsigned __int128 A0B0, A1B0_A0B1, tmplo;
	__int128 A1B1;
	
	A1B1 = (__int128) Ahi * Bhi;
	A0B0 = (__int128) Alo * Blo;
	A1B0_A0B1 = (__int128) ((__int128) ((__int128) Ahi * Blo) + ((__int128) Alo * Bhi)) + HIGH(A0B0);
	
	tmplo = (__int128) LOW(A0B0) | ((__int128)A1B0_A0B1 << 64);
	*Rhi += (__int128) A1B1 + HIGH(A1B0_A0B1) + add_overflow(Rlo, tmplo);
}

/*static inline void multadd128_karatsuba(__int128* restrict Rhi,*/
/*	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,*/
/*	const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R using karatsuba.*/
/*	unsigned __int128 A0B0, A1B0_A0B1;*/
/*	__int128 A1B1;*/
/*	*/
/*	A1B1 = (__int128) Ahi * Bhi;*/
/*	A0B0 = (__int128) Alo * Blo;*/
/*	A1B0_A0B1 = (__int128) ((__int128) ((__int128) A1B1 + A0B0) - ((__int128) ((__int128) Alo - Ahi) * ((__int128) Blo - Bhi))) + HIGH(A0B0);*/
/*	*/
/*	*Rhi += (__int128) A1B1 + HIGH(A1B0_A0B1) + add_overflow(Rlo, (__int128) LOW(A0B0) + ((__int128)A1B0_A0B1 << 64));*/
/*}*/

/*static inline void multadd128a(__int128* restrict Rhi,*/
/*	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,*/
/*	const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R, using mulx64.*/
/*	unsigned __int128 A0B0, tmplo;*/
/*	__int128 A1B1, A1B0, A0B1, aux1, aux2, aux3;*/
/*	*/
/*	uint64_t *hi, *lo;*/
/*	*/
/*	lo = (uint64_t*) &A1B1;*/
/*	hi = lo + 1;*/
/*	*lo = mulx64(Ahi, Bhi, hi);*/
/*	lo = (uint64_t*) &A1B0;*/
/*	hi = lo + 1;*/
/*	*lo = mulx64(Ahi, Blo, hi);*/
/*	lo = (uint64_t*) &A0B1;*/
/*	hi = lo + 1;*/
/*	*lo = mulx64(Alo, Bhi, hi);*/
/*	lo = (uint64_t*) &A0B0;*/
/*	hi = lo + 1;*/
/*	*lo = mulx64(Alo, Blo, hi);*/
/*	*/
/*	aux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);*/
/*	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);*/
/*	aux1 = (__int128) HIGH(A1B1);*/
/*	*/
/*	tmplo = (__int128) LOW(A0B0) | (aux3 << 64);*/
/*	*Rhi += (__int128) aux2 + (aux1 << 64) + __builtin_add_overflow(*Rlo, tmplo, Rlo);*/
/*}*/

/*static inline void multadd128g(__int128* restrict Rhi,*/
/*	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,*/
/*	const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R, using gmp's mpn functions.*/
/*	const uint64_t A[2] = { Alo, Ahi }, B[2] = { Blo, Bhi };*/
/*	uint64_t C[4] = {0}, R[4] = { LOW(*Rlo), HIGH(*Rlo), LOW(Rhi), HIGH(*Rhi)};*/
/*	//unsigned __int128 tmplo;*/
/*	*/
/*	*/
/*	mpn_mul_n(C, A, B, 2);*/
/*	*/
/*	mpn_add_n(R, R, C, 4);*/
/*	*/
/*	*Rlo = (__int128) R[0] | ((__int128) R[1] << 64);*/
/*	*Rhi = (__int128) R[2] | ((__int128) R[3] << 64);*/
/*}*/

/*void multadd128kara(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi,*/
/*	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R using karatsuba;*/
/*	unsigned __int128 A0B0, tmplo;*/
/*	__int128 A1B1, A1B0_A0B1, aux1, aux2, aux3;*/
/*	*/
/*	A1B1 = (__int128) Ahi * Bhi;*/
/*	A0B0 = (__int128) Alo * Blo;*/
/*	tmplo = (__int128) (Alo - Ahi) * (Blo - Bhi);*/
/*	//A1B0_A0B1 = (__int128) A0B0 + A1B1 - tmplo;*/
/*	aux1 = __builtin_add_overflow(A0B0, A1B1, &A1B0_A0B1);*/
/*	aux1 -= __builtin_sub_overflow(A1B0_A0B1, tmplo, &A1B0_A0B1);*/
/*	*/
/*	aux3 = (__int128) HI(A0B0) + LOW(A1B0_A0B1);*/
/*	aux2 = (__int128) HIGH(aux3) + HIGH(A1B0_A0B1) + LOW(A1B1);*/
/*	aux1 += (__int128) HIGH(A1B1);*/
/*	*/
/*	tmplo = (__int128) LOW(A0B0) | (aux3 << 64);*/
/*	*Rhi += (__int128) aux2 + (aux1 << 64) + add_overflow(Rlo, tmplo);*/
/*}*/

/*void m_multadd128(__int128* Rhi, unsigned __int128* Rlo, const uint64_t Ahi,*/
/*	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R for mult by M use;*/
/*	unsigned __int128 A0B0, tmplo;*/
/*	__int128 A1B1, A1B0_A0B1, aux1, aux2, aux3;*/
/*	*/
/*	A1B1 = (__int128) Ahi * Bhi;*/
/*	A0B0 = (__int128) Alo * Blo;*/
/*	aux3 = (__int128) (Alo - Ahi) * (Blo - Bhi);*/
/*	A1B0_A0B1 = (__int128) A0B0 + A1B1 - aux3;*/
/*	*/
/*	aux3 = (__int128) HI(A0B0) + LOW(A1B0_A0B1);*/
/*	aux2 = (__int128) HI(aux3) + HI(A1B0_A0B1) + LOW(A1B1);*/
/*	aux1 = (__int128) HI(A1B1);*/
/*	*/
/*	tmplo = *Rlo;*/
/*	*Rlo += (__int128) LOW(A0B0) + ((__int128)(HI(A0B0) + LOW(A1B0_A0B1)) << 64);*/
/*	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);*/
/*}*/

static inline void mm1_multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const uint64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R for mult by M or M1 use;
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0);
	
	tmplo = (__int128) A0B0 + ((__int128)(LOW(A0B1) + LOW(A1B0)) << 64);
	*Rhi += (__int128) aux2 + A1B1 + add_overflow(Rlo, tmplo);
}

/*static inline void mm1_multadd128__(__int128* restrict Rhi,*/
/*	unsigned __int128* restrict Rlo, const uint64_t Ahi, const uint64_t Alo,*/
/*	const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	*Rhi += (__int128) __builtin_add_overflow((__int128) Alo * Blo, *Rlo, Rlo)*/
/*		+ __builtin_add_overflow(*Rlo, ((unsigned __int128) ((__int128) Alo * Bhi) << 64), Rlo)*/
/*		+ __builtin_add_overflow(*Rlo, ((unsigned __int128) ((__int128) Ahi * Blo) << 64), Rlo)*/
/*		+ HIGH((__int128) Alo * Bhi) + HI((unsigned __int128) Ahi * Blo) + ((__int128) Ahi * Bhi);*/
/*}*/

/*static inline void mm1_multadd128x(__int128* restrict Rhi,*/
/*	unsigned __int128* restrict Rlo, const uint64_t Ahi, const uint64_t Alo,*/
/*	const int64_t Bhi, const uint64_t Blo)*/
/*{*/
/*	// multiplies A and B and adds the result to R for mult by M or M1 use;*/
/*	unsigned __int128 A0B0, A1B0, tmplo;*/
/*	__int128 aux1, aux2, aux3;*/
/*	*/
/*	A1B0 = (__int128) Ahi * Blo;*/
/*	A0B0 = (__int128) Alo * Blo;*/
/*	*/
/*	aux3 = (__int128) HIGH(A0B0) + LOW(Alo * Bhi) + LOW(Ahi * Blo);*/
/*	aux2 = (__int128) HIGH(aux3) + HIGH((__int128) Alo * Bhi) + HIGH(A1B0) + LOW(Ahi * Bhi);*/
/*	aux1 = (__int128) HIGH((__int128) Ahi * Bhi);*/
/*	*/
/*	tmplo = (__int128) LOW(Alo * Blo) | (aux3 << 64);*/
/*	*Rhi += (__int128) aux2 + (aux1 << 64) + __builtin_add_overflow(*Rlo, tmplo, Rlo);*/
/*}*/

static inline void OLD_m1_multadd128(unsigned __int128* restrict Rlo,
	const uint64_t Ahi, const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R for mult by M1 use;
	*Rlo += (__int128) Alo * Blo + ((__int128) (LOW(Alo * Bhi) + LOW(Ahi * Blo)) << 64);
}

static inline void m1_multadd128(unsigned __int128* restrict Rlo,
	unsigned __int128 A, __int128 B)
{
	// multiplies A and B and adds the result to R for mult by M1 use;
	*Rlo += (__int128) A * B;
}

inline void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	unsigned __int128 aux1, aux2;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		aux1 = (unsigned __int128) LOW(Rlo[i]) * LAMBDA;
		aux2 = (unsigned __int128) HI(Rlo[i]) * LAMBDA + HIGH(aux1);
		Rlo[i] = ((__int128) aux2 << 64) | aux1;
		Rhi[i] = (__int128) Rhi[i] * LAMBDA + HIGH(aux2);
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
		
		// Small trick for later
		Rhi[i] += Rlo[i] != 0;
	}
}

static inline void m_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		//Rhi[i] = 0;
		Rlo[i] = 0;
		for(j = 1; j < N - i; j++)
			mm1_multadd128(Rhi + i, Rlo + i, HIGH(A[i + j]), LOW(A[i + j]),
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			mm1_multadd128(Rhi + i, Rlo + i, HIGH(A[j]), LOW(A[j]),
				Mhi[i - j], Mlo[i - j]);
	}
}

static inline void m1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand (in pmns128 we only care about the lower 128 bits for
	// this operation).
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		for(j = 1; j < N - i; j++)
			m1_multadd128(Rlo + i, A[i+j], M1Lambda[N - j]);
				//((__int128) M1Lambdahi[N - j] << 64) | M1Lambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			m1_multadd128(Rlo + i, A[j], M1[i - j]);
				//((__int128) M1hi[i - j] << 64) | M1lo[i - j]);
	}
}

static inline void mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho.
	unsigned __int128 V[N];
	register uint16_t i;
	
	m1_mns128_mod_mult_ext_red(V, Rlo);
	
	m_mns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

inline void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	mns128_montg_int_red(res, Rhi, Rlo);
}

static inline void OLD_mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	unsigned __int128 aux;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		aux = (unsigned __int128) LOW(Rlo[i]) * (LAMBDA);
		aux = (unsigned __int128) HI(Rlo[i]) * (LAMBDA) + HIGH(aux);
		Rlo[i] = (__int128) Rlo[i] * (LAMBDA);
		Rhi[i] = (__int128) Rhi[i] * (LAMBDA) + HIGH(aux);
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

static inline void OLD_m_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, const restrict poly128 A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			mm1_multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			mm1_multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				Mhi[i - j], Mlo[i - j]);
	}
}

static inline void OLD_m1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	const restrict poly128 A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand (in pmns128 we only care about the lower 128 bits for
	// this operation).
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			OLD_m1_multadd128(Rlo + i, A->hi[i + j], A->lo[i + j],
				M1Lambdahi[N - j], M1Lambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			OLD_m1_multadd128(Rlo + i, A->hi[j], A->lo[j],
				M1hi[i - j], M1lo[i - j]);
	}
}

static inline void OLD_mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho.
	unsigned __int128 V[N], V2[N];
	register uint16_t i;
	
	
	memcpy(V, Rlo, N * sizeof(__int128));
	memcpy(V2, Rhi, N * sizeof(__int128));
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HI(Rlo[i]);
		Rlo[i] = 0;
	}
	
	OLD_m1_mns128_mod_mult_ext_red(Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HIGH(Rlo[i]);
		Rlo[i] = 0;
		Rhi[i] = 0;
	}
	
	OLD_m_mns128_mod_mult_ext_red(Rhi, Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		Rhi[i] = (__int128) V2[i] + Rhi[i] + (V[i] != 0);
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

void convert_string_to_amns128(restrict poly128 res, const char* string)
{
	uint8_t counter;
	register uint16_t i, j;
	const unsigned __int128 rho = ((__int128)1) << RHO;
	unsigned __int128 limb, Rlo[N] = {0};
	__int128 Rhi[N] = {0}, tmp[N] = {0};
	poly stok;
	init_poly(2 * N, &stok);
	
	convert_string_to_poly(&stok, string);
	
	if(stok->deg > 2 * N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	counter = 0;
	for(i = 0; i < N - 1; i++)
	{
		limb = (((unsigned __int128) stok->t[2 * i + 1]) << 64)
		 | ((uint64_t) stok->t[2 * i]);
		tmp[i] |= (limb << counter) & (rho - 1);
		tmp[i + 1] |= (limb >> (RHO - counter));
		counter = (counter + 128 - RHO) % RHO;
	}
	limb = (unsigned __int128) ((unsigned __int128) stok->t[2 * i]) |
			(((unsigned __int128) stok->t[2 * i + 1]) << 64);
	tmp[i] |= (limb << counter) & (rho - 1);
	
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			multadd128(Rhi + j, Rlo + j, HIGH(tmp[i]), LOW(tmp[i]),
				__Pihi__[i][j], __Pilo__[i][j]);
	
	mns128_montg_int_red(res, Rhi, Rlo);
	
end:
	free_poly(stok);
}

void convert_amns128_to_poly(restrict poly* res, const restrict poly128 P)
{
	register uint16_t i;
	poly aux, ag, tmp;
	poly128 a;
	__int128 Quitehi[N];
	unsigned __int128 Quitelo[N];
	
	init_polys(2 * N, &ag, &tmp, NULL);
	init_poly(2, &aux);
	init_poly128(N, &a);
	if((*res)->deg < 2 * N)
	{
		free_poly(*res);
		init_poly(2 * N, res);
	}
	(*res)->deg = 2 * N;
	for(i = 1; i < 2 * N; i++)
		(*res)->t[i] = 0;
	
	for(i = 0; i < N; i++)
	{
		Quitehi[i] = (__int128) P->hi[i];
		Quitelo[i] = (__int128) P->lo[i];
	}
	
	// THIS DOESN'T WORK
	// For an unknown reason, the result isn't times phi^-1 which disrupts things
	// Instead we use a workaround by multiplying by "one".
	// TODO: fix this problem
	mns128_montg_int_red(a, Quitehi, Quitelo);
	
	
	///////////////////////////// WORKAROUND ///////////////////////////////
	
/*	poly128 one;*/
/*	init_poly128(N, &one);*/
/*	for(i = 1; i < N; i++)*/
/*	{*/
/*		one->lo[i] = 0;*/
/*		one->hi[i] = 0;*/
/*	}*/
/*	one->hi[0] = 0;*/
/*	one->lo[0] = 1;*/
/*	*/
/*	amns128_montg_mult(a, P, one);*/
/*	*/
/*	//p128_print(a);*/
/*	*/
/*	free_poly128(one);*/
	
	///////////////////////////// WORKAROUND ///////////////////////////////
	
	(*res)->t[0] = a->lo[0];
	(*res)->t[1] = a->hi[0];
	for(i = 1; i < N; i++)
	{
		aux->t[0] = a->lo[i];
		aux->t[1] = a->hi[i];
		mp_mult(&ag, aux, &Gi[i - 1]);
		
		mp_copy(&tmp, *res);
		mp_add(res, tmp, ag);
	}
	
	mp_copy(&tmp, *res);
	mp_mod(res, tmp, &__P__);
	
	free_polys(aux, ag, tmp, NULL);
	free_poly128(a);
}


static inline void mns128_montg_int_red_pre(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho.
	unsigned __int128 V[N], V2[N];
	register uint16_t i;
	
	
	memcpy(V, Rlo, N * sizeof(__int128));
	memcpy(V2, Rhi, N * sizeof(__int128));
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HI(Rlo[i]);
		Rlo[i] = 0;
	}
	
	m1_mns128_mod_mult_ext_red_pre(Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HIGH(Rlo[i]);
		Rlo[i] = 0;
		Rhi[i] = 0;
	}
	
	m_mns128_mod_mult_ext_red_pre(Rhi, Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		Rhi[i] = (__int128) V2[i] + Rhi[i] + (V[i] != 0);
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

inline void amns128_montg_mult_pre(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red_pre(Rhi, Rlo, A, B);
	
	mns128_montg_int_red_pre(res, Rhi, Rlo);
}

inline void amns128_montg_mult_hyb(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	mns128_montg_int_red_pre(res, Rhi, Rlo);
}

static inline int64_t randomint64(void)
{
	// Function to generate a random 64 bit number.
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrhohi(int64_t param)
{
	// Utility function to get a high part for a random poly128.
	return param & ((1ULL<<(RHO - 64)) - 1);
}

void randpoly128(poly128 P)
{
	// Generates a random poly128 with appropriate, lower than rho high part.
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->lo[i] = randomint64();
		if(RHO > 64) P->hi[i] = __modrhohi(randomint64()) * (1 + (rand() & 1) * -2);
	}
}

void __main__(void)
{
/*	int64_t Ahi = 0x5c096e6b558ba549, Bhi = 0x5918460d4a05af9c;*/
/*	uint64_t  Alo = 0xf9dadbd740482391, Blo = 0xc4703ac88b18e4c3;*/
/*	*/
/*	const char RES[] = */
/*	"0x2008017506164ba2bbc36636e65b8f4caca8fc6b0e76b6a3dd70ff9147383b73";*/
	
/*	int64_t Ahi = 0xfc096e6b558ba549, Bhi = 0xf918460d4a05af9c;*/
/*	uint64_t  Alo = 0xf9dadbd740482391, Blo = 0xc4703ac88b18e4c3;*/
/*	*/
/*	const char RES[] = */
/*	"0x1b5dc7ca3fcbcc34673dbafa172c2d2ca8fc6b0e76b6a3dd70ff9147383b73";*/
	
/*	int64_t Ahi = 0xffffffffffffffff, Bhi = 0xffffffffffffffff;*/
/*	uint64_t  Alo = 0xffffffffffffffff, Blo = 0xffffffffffffffff;*/
/*	*/
/*	const char RES[] = */
/*	"0x0000000000000000000000000000000000000000000000000";*/
	
/*	int64_t Ahi = 0x7e4deecd7fb21228, Bhi = -0x11fb30a7f4f269dd;*/
/*	uint64_t  Alo = 0xbe23c803cef47a1, Blo = 0xc97d66ba9ffee37;*/
/*	*/
/*	const char RES[] = */
/*	"0xf720e4b9bc91c7b14365a9592434dee8edbb60fd615d17fd17ab67202e5f1197";*/
	
/*	int64_t Ahi = 0xead7e6d770b694af, Bhi = 0xa3c338cf5fdb7ee7;*/
/*	uint64_t  Alo = 0xcd3ad86ad0c870de, Blo = 0x83ed23c6860ab850;*/
	
	const char RES[] =
		"0x23fbe09528a6b4ada3fa9544c0b0b34f51d40d5578d6c57b14af738e2a3c5388";
	
/*	__int128 R1 = 0;*/
/*	unsigned __int128 R2 = 0;*/
	
	__int128 R1 = (__int128) ((__int128) 0x8dc15207aec5564b << 64) + 0x4ba591887a57f6ed;
	unsigned __int128 R2 = (__int128) ((__int128) 0x4b5f46c669762ca7 << 64) + 0x4f4e53803dcd7e28;
	
/*	__int128 R1 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;*/
/*	unsigned __int128 R2 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;*/
	
	//multadd128g(&R1, &R2, Ahi, Alo, Bhi, Blo);
	
	printf("%s\n", RES);
	
	printf("0x%lx%016lx%016lx%016lx\n", HIGH(R1), LOW(R1), HIGH(R2), LOW(R2));
	
/*	printf("\n\n");*/
	
/*	poly128 A;*/
/*	init_poly128(N, &A);*/
/*	const char a[] = "bf6dc9f34905d4ccea18b34313d7f22412795efa0161f7ebcd5912a900ea7d255661bb894729e4fd85a477d3c575f3e97fcd1e6e2fd01d5317724f38def3c7f944162bb4ae4dcd5b1522efca1f3713a927c91f1113096ced7585edf7fef8cc9334dc56e8483a3c49f4a0fb9bb73c00b8b00e3d11435184eacbd45dd38fcbcadd";*/
/*	*/
/*	convert_string_to_amns128(A, a);*/
/*	*/
/*	p128_print(A);*/
/*	*/
/*	free_poly128(A);*/
}

void __benchmult__(void)
{
	uint64_t sum1, sum2, c = clock(), lo1, lo2;
	int64_t hi1, hi2;
	__int128 dummy11 = 0, dummy21 = 0;
	unsigned __int128 dummy12 = 0, dummy22 = 0;
	
	srand((unsigned) (time(&hi1)));
	
	sum1 = 0;
	sum2 = 0;
	for(int i = 0; i < 100000; i++)
	{
		hi1 = randomint64();
		hi2 = randomint64();
		lo1 = randomint64();
		lo2 = randomint64();
		c = clock();
		multadd128x(&dummy11, &dummy12, hi1, lo1, hi2, lo2);
		sum1 += clock() - c;
		c = clock();
		multadd128(&dummy21, &dummy22, hi1, lo1, hi2, lo2);
		sum2 += clock() - c;
	}
	sum1 = 0;
	sum2 = 0;
	for(int i = 0; i < 100000; i++)
	{
		hi1 = randomint64();
		hi2 = randomint64();
		lo1 = randomint64();
		lo2 = randomint64();
		c = clock();
		multadd128x(&dummy11, &dummy12, hi1, lo1, hi2, lo2);
		sum1 += clock() - c;
		c = clock();
		multadd128(&dummy21, &dummy22, hi1, lo1, hi2, lo2);
		sum2 += clock() - c;
	}
	printf("1: %ld\n2: %ld\n", sum1, sum2);
	__print128(dummy11);
	__print128(dummy21);
	__print128(dummy12);
	__print128(dummy22);
}

void __multchecks__(void)
{
	poly128 a, b, c;
	init_poly128s(N, &a, &b, &c, NULL);
	int64_t seed;
	
	srand((unsigned) (time(&seed)));
	
	for(int i = 0; i < 100; i++)
	{
		randpoly128(a);
		randpoly128(b);
		p128_print(a);
		p128_print(b);
		amns128_montg_mult(c, a, b);
		p128_print(c);
	}
	
	free_poly128s(a, b, c, NULL);
}

void __multbench__(void)
{
	uint64_t c1 = 0, sum;
	poly128 a, b, c, soak1, soak2;
	init_poly128s(N, &a, &b, &c, &soak1, &soak2, NULL);
	
	srand((unsigned) (time(((int64_t*)(&c1)))));
	
	c1 = clock();
	
	randpoly128(soak2);
	soak2->lo[0] += Gi[0].t[0];
	soak2->lo[0] += __P__.t[0];
	
	for(int i = 0; i < 1000; i++)
	{
		randpoly128(a);
		randpoly128(b);
		amns128_montg_mult(c, a, b);
		amns128_montg_mult(soak1, c, soak2);
		amns128_montg_mult(soak2, c, soak1);
	}
	
	sum = 0;
	for(int i = 0; i < 10000; i++)
	{
		randpoly128(a);
		randpoly128(b);
		c1 = clock();
		amns128_montg_mult(c, a, b);
		sum += clock() - c1;
		amns128_montg_mult(soak1, c, soak2);
		amns128_montg_mult(soak2, c, soak1);
	}
	
	printf("%ld\n", sum);
	
	free_poly128s(a, b, c, soak1, soak2, NULL);
}

void __full_mult128_demo(void)
{
	const char a[] = "eedad1f0fcbd8935abfccb2a2752c672db377b0575464c626bd4003dca25c7c80589ac8e48c02c2379c73d45d9558ba5e708e338ec0686a4ab1356e9ad34348859398e506895f6c321bc4f8d2e2f21783ff2e65030a968331c9829e17900a1af32f1b7a3131ab61c0ad9b9e98a23f932d2f8936bc32d06128b086e619a1e35c1",
		b[] = "a3442cfd10dcb1d9a20a8541a8551a37d23bbe02dff3faeaf4fd5253975bfad4c759b9101bf4ae5197a18385f0658f5e436644638cd923ef490a7251402fdc0ef87dd61af1f0dbc9b81ef7de744b062a0d5df758eb1c1c41fd4cf01255aba467f1f673c4e5b20540f9c77ae0e430a323225346fdcc804912a09c6728b5837935",
		c[] = "f85d170fe1a19c62f1c527086d1166ad84330e32b2845c04c3241b5bd7ff0806da511a61e67ec3a8ab988f035ad03437ed79f85ea02d657f73994e62ca71b42ef0fbeb8cdc1ac744e96110dcb69af62a34b830814845938870e396e641ce6729ec6953e6639cc7ff0bd4fe4d34452caa974673fe52006821764b8e456d51779a";
	
	poly128 A, B, C;
	poly aux;
	
	init_poly128s(N, &A, &B, &C, NULL);
	init_poly(2 * N, &aux);
	
	
	convert_string_to_amns128(A, a);
	convert_string_to_amns128(B, b);
	
/*	printf("0x%s\n", a);*/
/*	p128_print(A);*/
/*	printf("\n0x%s\n", b);*/
/*	p128_print(B);*/
	
	convert_string_to_amns128(C, "1");
	p128_print(C);
	p128_print(&__theta__);
	
/*	amns128_montg_mult(C, A, B);*/
	
/*	p128_print(C);*/
	
/*	convert_amns128_to_poly(&aux, C);*/
	
/*	printf("\n%s\n", c);*/
/*	mp_print(aux);*/
	
	free_poly128s(A, B, C, NULL);
	free_poly(aux);
}

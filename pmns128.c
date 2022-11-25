#include <stdio.h>
#include <time.h>
//#include <gmp.h>
#include <string.h>

#include "pmns128.h"
#include "utilitymp.h"


static inline void multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const uint64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// Multiplies A and B and adds the result to R using the classic
	// schoolbook algorithm
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

static inline void mns_multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// Multiplies A and B and adds the result to R faster in the case where A
	// and B are both in the PMNS and therefore their coefficients are lower
	// than rho. In that case the middle sum will never overflow, speeding
	// up the calculations.
	unsigned __int128 A0B0, A1B0_A0B1, tmplo;
	__int128 A1B1;
	
	A1B1 = (__int128) Ahi * Bhi;
	A0B0 = (__int128) Alo * Blo;
	A1B0_A0B1 = (__int128) ((__int128) ((__int128) Ahi * Blo) + ((__int128) Alo * Bhi)) + HIGH(A0B0);
	
	tmplo = (__int128) LOW(A0B0) | ((__int128)A1B0_A0B1 << 64);
	*Rhi += (__int128) A1B1 + HIGH(A1B0_A0B1) + add_overflow(Rlo, tmplo);
}

static inline void m1_multadd128(unsigned __int128* restrict Rlo,
	unsigned __int128 A, __int128 B)
{
	// Multiplies A and B and adds the result to R in the case where we only
	// want the lower part of the result.
	*Rlo += (__int128) A * B;
}

#ifdef LAMBDA

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
			mns_multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		aux1 = (unsigned __int128) LOW(Rlo[i]) * LAMBDA;
		aux2 = (unsigned __int128) HI(Rlo[i]) * LAMBDA + HIGH(aux1);
		Rlo[i] = ((__int128) aux2 << 64) | LOW(aux1);
		Rhi[i] = (__int128) Rhi[i] * LAMBDA + HIGH(aux2);
		
		for(j = 0; j < i + 1; j++)
			mns_multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

#endif

#ifdef LENEXTPOLY

inline void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) an irreducible polynomial used for reduction. Result in R
	register uint16_t i, j, k;
	__int128 Thi;
	unsigned __int128 aux1, aux2, Tlo;
	
	for(i = 0; i < N; i++)
	{
		Thi = 0;
		Tlo = 0;
		for(j = 1; j < N - i; j++)
			mns_multadd128(&Thi, &Tlo, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		for(k = 0; (k < LENEXTPOLY) && (i + k < N); k++)
		{
			aux1 = (unsigned __int128) LOW(Tlo) * EXTPOLY[k] + LOW(Rlo[i + k]);
			aux2 = (unsigned __int128) HI(Tlo) * EXTPOLY[k] + HIGH(aux1) + HI(Rlo[i + k]);
			Rlo[i + k] = ((__int128) aux2 << 64) | LOW(aux1);
			Rhi[i + k] += (__int128) Thi * EXTPOLY[k] + HIGH(aux2);
		}
		
		for(j = 0; j < i + 1; j++)
			mns_multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

#endif

#ifdef M_or_B_is_M

static inline void m_or_b_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		Rhi[i] += Rlo[i] != 0;
		Rlo[i] = 0;
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, HIGH(A[i + j]), LOW(A[i + j]),
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, HIGH(A[j]), LOW(A[j]),
				Mhi[i - j], Mlo[i - j]);
	}
}

static inline void m1_or_b1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R. 
	// In pmns128 we only care about the lower 128 bits for this operation.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		for(j = 1; j < N - i; j++)
			m1_multadd128(Rlo + i, A[i+j], M1Lambda[N - j]);
		
		for(j = 0; j < i + 1; j++)
			m1_multadd128(Rlo + i, A[j], M1[i - j]);
	}
}

#endif

#ifdef M_or_B_is_B

static inline void m_or_b_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		Rhi[i] += Rlo[i] != 0;
		Rlo[i] = 0;
		for(j = 0; j < N; j++)
			multadd128(Rhi + i, Rlo + i, HIGH(A[j]), LOW(A[j]),
				Bhi[j][i], Blo[j][i]);
	}
}

static inline void m1_or_b1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Vector-Matrix multiplication between A and B1, result in R.
	// In 128 bit version we only care about the lower 128 bits.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		for(j = 0; j < N; j++)
			m1_multadd128(Rlo + i, A[j], B1[j][i]);
	}
}

#endif

static inline void mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho using Montgomery's internal reduction algorithm.
	unsigned __int128 V[N];
	register uint16_t i;
	
	m1_or_b1_mns128_mod_mult_ext_red(V, Rlo);
	
	m_or_b_mns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

static inline void UNROLLED_mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Unrolled version.
	unsigned __int128 V[N];
	register uint16_t i;
	
	UNROLLED_m1_or_b1_mns128_mod_mult_ext_red(V, Rlo);
	
	UNROLLED_m_or_b_mns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

inline void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	mns128_montg_int_red(res, Rhi, Rlo);
}

inline void UNROLLED_amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Unrolled version.
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	UNROLLED_mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	UNROLLED_mns128_montg_int_red(res, Rhi, Rlo);
}

void amns128_rtl_sqandmult(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Right to Left version. Slower than ltr.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly128 tmp;
	init_poly128(N, &tmp);
	
	p128_copy(res, &__theta__);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		for(j = 0; j < 64; j++)
			{
				if(aux & (1ULL << j))
					amns128_montg_mult(tmp, res, base);
				else
					amns128_montg_mult(tmp, res, &__theta__);
				amns128_montg_mult(res, tmp, tmp);
			}
	}
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			amns128_montg_mult(tmp, res, base);
		else
			amns128_montg_mult(tmp, res, &__theta__);
		amns128_montg_mult(res, tmp, tmp);
		aux >>= 1;
	}
	
	free_poly128(tmp);
}

void amns128_ltr_sqandmult(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Left to Right version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly128 tmp;
	init_poly128(N, &tmp);
	
	p128_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		amns128_montg_mult(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			amns128_montg_mult(res, tmp, base);
		else
			amns128_montg_mult(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			amns128_montg_mult(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				amns128_montg_mult(res, tmp, base);
			else
				amns128_montg_mult(res, tmp, &__theta__);
		}
	}
	
	free_poly128(tmp);
}

void amns128_montg_ladder(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux, b;
	
	poly128 tmp, R[2];
	init_poly128(N, &tmp);
	R[0] = res;
	R[1] = tmp;
	
	p128_copy(res, &__theta__);
	p128_copy(tmp, base);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		amns128_montg_mult(R[1 - b], R[1 - b], R[b]);
		amns128_montg_mult(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			amns128_montg_mult(R[1 - b], R[1 - b], R[b]);
			amns128_montg_mult(R[b], R[b], R[b]);
		}
	}
	
	free_poly128(tmp);
}

inline void amns128_sqandmult(restrict poly128 res, const char* restrict base,
	const char* restrict exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Uses montgomery ladder.
	
	poly128 pbase;
	mpnum mpexponent;
	init_poly128(N, &pbase);
	init_mpnum(1, &mpexponent);
	
	convert_string_to_amns128(pbase, base);
	convert_string_to_multipre(&mpexponent, exponent);
	
	amns128_montg_ladder(res, pbase, mpexponent);
	
	free_poly128(pbase);
	free_mpnum(mpexponent);
}

void convert_string_to_amns128(restrict poly128 res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (we multiply it by PHI in the
	// process).
	
	uint8_t counter;
	register uint16_t i, j;
	const unsigned __int128 rho = ((__int128)1) << RHO;
	unsigned __int128 limb, Rlo[N] = {0};
	__int128 Rhi[N] = {0}, tmp[N] = {0};
	mpnum stok, aux1, aux2;
	init_mpnums(2 * N, &stok, &aux1, &aux2, NULL);
	
	convert_string_to_multipre(&aux1, string);
	
	mp_mod(&aux2, aux1, &__P__);
	
	for(i = 0; i < aux2->deg; i++)
		stok->t[i] = aux2->t[i];
	
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
	
	free_mpnums(stok, aux1, aux2, NULL);
}

void convert_amns128_to_multipre(restrict mpnum* res, const restrict poly128 P)
{
	// Function that converts out of the AMNS system and into a multiprecision
	// number.
	
	register uint16_t i;
	mpnum aux, ag, tmp, tmp2;
	poly128 a;
	
	init_mpnums(2 * N, &ag, &tmp, &tmp2, NULL);
	init_mpnum(2, &aux);
	init_poly128(N, &a);
	if((*res)->deg < 2 * N)
	{
		free_mpnum(*res);
		init_mpnum(2 * N, res);
	}
	(*res)->deg = 2 * N;
	for(i = 1; i < 2 * N; i++)
		(*res)->t[i] = 0;
	
	poly128 one;
	init_poly128(N, &one);
	for(i = 1; i < N; i++)
	{
		one->lo[i] = 0;
		one->hi[i] = 0;
	}
	one->hi[0] = 0;
	one->lo[0] = 1;
	
	amns128_montg_mult(a, P, one);
	
	free_poly128(one);
	
	(*res)->sign = 1 - 2 * (a->hi[0] < 0);
	(*res)->t[0] = a->lo[0] * (*res)->sign;
	(*res)->t[1] = a->hi[0] * (*res)->sign - ((*res)->sign == -1);
	for(i = 1; i < N; i++)
	{
		aux->sign = 1 - 2 * (a->hi[i] < 0);
		aux->t[0] = a->lo[i] * aux->sign;
		aux->t[1] = a->hi[i] * aux->sign - (aux->sign == -1);
		
		mp_mult(&ag, aux, &Gi[i - 1]);
		
		mp_copy(&tmp, *res);
		mp_add(res, tmp, ag);
	}
	mp_copy(&tmp, *res);
	mp_mod(res, tmp, &__P__);
	
	free_mpnums(aux, ag, tmp, tmp2, NULL);
	free_poly128(a);
}

static inline int64_t randomint64(void)
{
	// Function to generate a random 64 bit number.
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrhohi(int64_t param)
{
	// Utility function to get a high part for a random poly128.
	if (RHO <= 64)
		return 0;
	else
		return param & ((1ULL << (RHO - 64)) - 1);
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

void __multchecks__(char* nbmults)
{
	// Debug tool to quickly check if the pmns gives the correct results.
	poly128 a, b, c;
	init_poly128s(N, &a, &b, &c, NULL);
	int64_t seed;
	
	register uint64_t cap = 100;
	
	if(nbmults[0] != '\0')
		cap = atoll(nbmults);
	
	srand((unsigned) (time(&seed)));
	
	for(register uint64_t i = 0; i < cap; i++)
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

void __sqandmultdemo(void)
{
	// TODO: delete it, used to check point to point process.
	const char a[] = "0xffffffffffffffffffffffffffffffffffffffffffffffff",
		b[] = "0xffffffffffffffffffffffffffffffffffffffffffffffff",
		c[] = "0xfca02b1ec94f1d1a5afdbacd6c0b04b9ded8463b62fade4d35e64f3de3a20e8b172382d318083548be85676b8b7dd347c32841c3efa88c928e7cc945f5e2253708053c8f75ff222717ea374b7f0b6ceafd6ac0b119d23300ac9b6a23cd0b332b438f7be9f6c4ee4cdbbe93444343a8011583e10a8c4cc0f07b67d1467b9a73e";
	
	poly128 C;
	mpnum aux;
	init_poly128(N, &C);
	init_mpnum(1, &aux);
	
	amns128_sqandmult(C, a, b);
	
	convert_amns128_to_multipre(&aux, C);
	
	printf("\n%s\n\n", c);
	mp_print(aux);
	printf("\n");
	
	free_poly128(C);
	free_mpnum(aux);
}

#include <stdio.h>
#include <time.h>
/*#include <omp.h>*/

#include "pmns.h"
#include "utilitymp.h"

#define NBTHREADZ 8

// -fopenmp
// _PRAGMAGCCUNROLLLOOP_ private(j, k)

#ifdef LAMBDA
inline void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	register uint16_t i, j;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
		
		R[i] = R[i] * LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

/*inline void karamns_mod_mult_ext_red(__int128* restrict R,*/
/*	const restrict poly A, const restrict poly B)*/
/*{*/
/*	// Function that multiplies A by B and applies external reduction using*/
/*	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.*/
/*	// Karatsuba version.*/
/*	*/
/*	register uint16_t i, j;*/
/*	*/
/*	__int128 Lo[N/2] = {0}, Hi[N] = {0}, Mid[N] = {0};*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*	{*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = 0; j < i + 1; j++)*/
/*		{*/
/*			Lo[i] += (__int128) A->t[i - j] * B->t[j];*/
/*			Mid[i] += (__int128) (A->t[i - j] + A->t[i - j + N/2]) * (B->t[j] + B->t[j + N/2]);*/
/*			Hi[i] += (__int128) A->t[N/2 + i - j] * B->t[N/2 + j];*/
/*		}*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = i + 1; j < N/2; j++)*/
/*		{*/
/*			Hi[i] -= (__int128) A->t[i + N/2 - j] * B->t[j];*/
/*			Hi[i + N/2] += (__int128) A->t[N + i - j] * B->t[N/2 + j];*/
/*			Mid[i + N/2] += (__int128) (A->t[i + N/2 - j] + A->t[i - j + N]) * (B->t[j] + B->t[j + N/2]);*/
/*		}*/
/*	}*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i + N/2] +=  Hi[i + N/2] * LAMBDA + Mid[i] - Lo[i] - Hi[i];*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i] += Lo[i] + (Mid[i + N/2] + Hi[i] - Hi[i + N/2]) * LAMBDA;*/
/*}*/

#ifdef M_or_B_is_M

static inline void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) (A[i + j]) * M[N - j] * LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) (A[j]) * M[i - j];
	}
}

/*static inline void karam_or_b_mns_mod_mult_ext_red(__int128* restrict R,*/
/*	int64_t* restrict A)*/
/*{*/
/*	// Function that multiplies A by M and applies external reduction using*/
/*	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.*/
/*	// Karatsuba version.*/
/*	*/
/*	register uint16_t i, j;*/
/*	*/
/*	__int128 Lo[N] = {0}, Hi[N] = {0}, Mid[N] = {0};*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*	{*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = 0; j < i + 1; j++)*/
/*		{*/
/*			Lo[i] += (__int128) A[i - j] * M[j];*/
/*			Hi[i] += (__int128) A[N/2 + i - j] * M[N/2 + j];*/
/*			Mid[i] += ((__int128)A[i - j] + A[i - j + N/2]) * (M[j] + M[j + N/2]);*/
/*		}*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = i + 1; j < N/2; j++)*/
/*		{*/
/*			Lo[i + N/2] += (__int128) A[i + N/2 - j] * M[j];*/
/*			Hi[i + N/2] += (__int128) A[N + i - j] * M[N/2 + j];*/
/*			Mid[i + N/2] += ((__int128)A[i + N/2 - j] + A[i - j + N]) * (M[j] + M[j + N/2]);*/
/*		}*/
/*	}*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i + N/2] +=  Lo[i + N/2] + Hi[i + N/2] * LAMBDA + Mid[i] -  Lo[i] -  Hi[i];*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i] += Lo[i] + (Mid[i + N/2] + Hi[i] - Lo[i + N/2] -  Hi[i + N/2]) * LAMBDA;*/
/*}*/



static inline void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			R[i] += ((uint64_t)A[i + j]) * M1[N - j] * LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += ((uint64_t)A[j]) * M1[i - j];
	}
}

/*static inline void karam1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,*/
/*	__int128* restrict A)*/
/*{*/
/*	// Function that multiplies A by M1 and applies external reduction using*/
/*	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.*/
/*	// Karatsuba version.*/
/*	*/
/*	register uint16_t i, j;*/
/*	*/
/*	uint64_t Lo[N] = {0}, Hi[N] = {0}, Mid[N] = {0};*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*	{*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = 0; j < i + 1; j++)*/
/*		{*/
/*			Lo[i] += A[i - j] * M1[j];*/
/*			Hi[i] += A[N/2 + i - j] * M1[N/2 + j];*/
/*			Mid[i] += (A[i - j] + A[i - j + N/2]) * (M1[j] + M1[j + N/2]);*/
/*		}*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = i + 1; j < N/2; j++)*/
/*		{*/
/*			Lo[i + N/2] += A[i + N/2 - j] * M1[j];*/
/*			Hi[i + N/2] += A[N + i - j] * M1[N/2 + j];*/
/*			Mid[i + N/2] += (A[i + N/2 - j] + A[i - j + N]) * (M1[j] + M1[j + N/2]);*/
/*		}*/
/*	}*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i + N/2] = Lo[i + N/2] + Hi[i + N/2] * LAMBDA + Mid[i] - Lo[i] - Hi[i];*/
/*	*/
/*	_PRAGMAGCCUNROLLLOOP_*/
/*	for(i = 0; i < N/2; i++)*/
/*		R[i] = Lo[i] + (Mid[i + N/2] + Hi[i] - Lo[i + N/2] -  Hi[i + N/2]) * LAMBDA;*/
/*}*/

#endif

#endif

#ifdef LENEXTPOLY

inline void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) an irreducible polynomial used for reduction. Result in R
	register uint16_t i, j, k;
	__int128 T;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			T += (__int128) A->t[i + j] * B->t[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

#ifdef M_or_B_is_M

static inline void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j, k;
	__int128 T;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			T += (__int128) A[i + j] * M[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * M[i - j];
	}
}

static inline void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j, k;
	__int128 T;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 1; j < N - i; j++)
			T += ((uint64_t)A[i + j]) * M1[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
			R[i] += ((uint64_t)A[j]) * M1[i - j];
	}
}

#endif

#endif


#ifdef M_or_B_is_B

static inline void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	register uint16_t i, j;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < N; j++)
			R[i] += (__int128) A[j] * B[j][i];
	}
}

static inline void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Vector-Matrix multiplication between A and B1, result in R.
	register uint16_t i, j;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < N; j++)
			R[i] += ((uint64_t)A[j]) * B1[j][i];
	}
}

#endif

inline void mns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Montgomery method.
	int64_t T[N];
	register uint16_t i;
	
	m1_or_b1_mns_mod_mult_ext_red(T, R);
	
	m_or_b_mns_mod_mult_ext_red(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

inline void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	mns_mod_mult_ext_red(R, A, B);

	mns_montg_int_red(res, R);
}

void amns_rtl_sqandmult(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Right to Left version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp, num, stok;
	init_polys(N, &tmp, &num, &stok, NULL);
	
	poly_copy(res, &__theta__);
	poly_copy(stok, base);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		// We do it in pairs to avoid the use of poly_copy
		for(j = 0; j < 32; j++)
			{
				if(aux & 1)
					amns_montg_mult(tmp, res, stok);
				else
					amns_montg_mult(tmp, res, &__theta__);
				amns_montg_mult(num, stok, stok);
				aux >>= 1;
				if(aux & 1)
					amns_montg_mult(res, tmp, num);
				else
					amns_montg_mult(res, tmp, &__theta__);
				amns_montg_mult(stok, num, num);
				aux >>= 1;
			}
	}
	// If our exponent has an odd number of bits we do this one time too many.
	// For this reason unless we find another way, Left to Right is preferred.
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			amns_montg_mult(tmp, res, stok);
		else
			amns_montg_mult(tmp, res, &__theta__);
		amns_montg_mult(num, stok, stok);
		aux >>= 1;
		if(aux & 1)
			amns_montg_mult(res, tmp, num);
		else
			amns_montg_mult(res, tmp, &__theta__);
		amns_montg_mult(stok, num, num);
		aux >>= 1;
	}
	
	free_polys(tmp, num, stok, NULL);
}

void amns_ltr_sqandmult(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Left to Right version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		amns_montg_mult(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			amns_montg_mult(res, tmp, base);
		else
			amns_montg_mult(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			amns_montg_mult(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				amns_montg_mult(res, tmp, base);
			else
				amns_montg_mult(res, tmp, &__theta__);
		}
	}
	
	free_poly(tmp);
}

void amns_montg_ladder(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux, b;
	
	poly tmp, R[2];
	init_poly(N, &tmp);
	R[0] = res;
	R[1] = tmp;
	
	poly_copy(res, &__theta__);
	poly_copy(tmp, base);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		amns_montg_mult(R[1 - b], R[1 - b], R[b]);
		amns_montg_mult(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			amns_montg_mult(R[1 - b], R[1 - b], R[b]);
			amns_montg_mult(R[b], R[b], R[b]);
		}
	}
	
	free_poly(tmp);
}

inline void amns_sqandmult(restrict poly res, const char* restrict base,
	const char* restrict exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Uses Left to Right square and multiply.
	
	poly pbase;
	mpnum mpexponent;
	init_poly(N, &pbase);
	init_mpnum(1, &mpexponent);
	
	convert_string_to_amns(pbase, base);
	convert_string_to_multipre(&mpexponent, exponent);
	
	amns_ltr_sqandmult(res, pbase, mpexponent);
	
	free_poly(pbase);
	free_mpnum(mpexponent);
}

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

void randpoly(poly P)
{
	// Function that generates a random polynomial within our PMNS system.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64()) * (1 + (rand() & 1) * -2);
}

void convert_string_to_amns(restrict poly res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (we multiply it by PHI in the
	// process).
	
	uint8_t counter;
	register uint16_t i, j;
	const uint64_t rho = (1ULL<<RHO);
	__int128 R[N] = {0};
	mpnum stok, tmp;
	init_mpnums(N, &stok, &tmp, NULL);
	
	convert_string_to_multipre(&tmp, string);
	
	mp_mod(&stok, tmp, &__P__);
	
	for(i = 1; i < N; i++)
		res->t[i] = 0;
	
	counter = 0;
	res->t[0] = stok->t[0] & (rho - 1);
	for(i = 1; i < stok->deg; i++)
	{
		counter = (counter + RHO) % 64; 
		res->t[i] = (stok->t[i-1] >> counter);
		res->t[i] |= (stok->t[i] << (64 - counter));
		res->t[i] &= (rho - 1);
	}
	counter = (counter + RHO) % 64;
	res->t[i] = (stok->t[i-1] >> counter);
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);
	
	free_mpnums(stok, tmp, NULL);
}

void convert_amns_to_multipre(restrict mpnum* res, const restrict poly P)
{
	// Function that converts out of the AMNS system and into a multiprecision
	// number.
	
	register uint16_t i;
	poly a;
	mpnum aux, ag, tmp, tmp2;
	__int128 Quite[N];
	
	init_poly(N, &a);
	init_mpnums(N, &ag, &tmp, &tmp2, NULL);
	init_mpnum(1, &aux);
	if((*res)->deg < N)
	{
		free_mpnum(*res);
		init_mpnum(N, res);
	}
	(*res)->deg = N;
	for(i = 1; i < N; i++)
	{
		(*res)->t[i] = 0;
		Quite[i] = (__int128) P->t[i];
	}
	Quite[0] = (__int128) P->t[0];
	
	mns_montg_int_red(a, Quite);
	
	// We use Horner to get out of the PMNS for now as the other method has issues.
	(*res)->sign = 1 - 2 * (a->t[N - 1] < 0);
	(*res)->t[0] = a->t[N - 1] * (*res)->sign;
	for(i = 0; i < N - 1; i++)
	{
		mp_mult(&tmp, *res, &Gi[0]);
		aux->sign = 1 - 2 * (a->t[N - 2 - i] < 0);
		aux->t[0] = a->t[N - 2 - i] * aux->sign;
		mp_add(&tmp2, tmp, aux);
		mp_mod(res, tmp2, &__P__);
	}
	
/*	(*res)->t[0] = a->t[0];*/
/*	for(i = 1; i < N; i++)*/
/*	{*/
/*		aux->sign = 1 - 2 * (a->t[i] < 0);*/
/*		aux->t[0] = a->t[i] * aux->sign;*/
/*		mp_mult(&ag, aux, &Gi[i - 1]);*/
/*		mp_copy(&tmp, *res);*/
/*		mp_add(res, tmp, ag);*/
/*	}*/
/*	*/
/*	mp_copy(&tmp, *res);*/
/*	mp_mod(res, tmp, &__P__);*/
	
	free_mpnums(aux, ag, tmp, tmp2, NULL);
	free_poly(a);
}

void __multchecks__(char* nbmults)
{
	// Used as a debug tool to see if the PMNS correctly gives us the proper
	// results with a few random values.
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	int64_t seed;
	
	register uint64_t cap = 100;
	
	if(nbmults[0] != '\0')
		cap = atoll(nbmults);
	
	srand((unsigned) (time(&seed)));
	
	for(register uint64_t i = 0; i < cap; i++)
	{
		randpoly(a);
		randpoly(b);
		poly_print(a);
		poly_print(b);
		amns_montg_mult(c, a, b);
		poly_print(c);
	}
	
	free_polys(a, b, c, NULL);
}

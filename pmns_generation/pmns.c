#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#include "structs.h"
#include "utilitymp.h"
#include "params.h"

#define SCHOOLBOOK(X) {\
	for(int i = 0; i < X; i++)\
	{\
		rop[i] = 0;\
		for(int j = 0; j < X; j++)\
			rop[i] += (__int128) vect[j] * matr[X - 1 - j + i];\
	}\
}

#define TOEP22TOP(X, F) {\
	__int128 t0[X/2], t1[X/2], t2[X/2];\
	int64_t v0p1[X/2], m0m1[X-1], m0m2[X-1];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = vect[i] + vect[i + X/2];\
	}\
	for(int i = 0; i < X-1; i++)\
	{\
		m0m1[i] = matr[i + X/2] - matr[i + X];\
		m0m2[i] = matr[i + X/2] - matr[i];\
	}\
	F (t0, v0p1, matr + X/2);\
	F (t1, vect, m0m1);\
	F (t2, vect + X/2, m0m2);\
	for(int i = 0; i < X/2; i++)\
	{\
		rop[i] = t0[i] - t2[i];\
		rop[i + X/2] = t0[i] - t1[i];\
	}\
}

#define TOEP33TOP(X, F) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3];\
	int64_t m03, v1m2[X/3], v0m2[X/3], v0m1[X/3], m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1], m0[(2*X/3) - 1], m1[(2*X/3) - 1], m3[(2*X/3) - 1];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = vect[i] - vect[X/3 + i];\
		v0m2[i] = vect[i] - vect[2*X/3 + i];\
	}\
	for(int i = 0; i < (2*X/3) - 1; i++)\
	{\
		m0[i] = matr[i + 2*X/3];\
		m1[i] = matr[i + X];\
		m3[i] = matr[i + X/3];\
		m03 = m0[i] + m3[i];\
		m034[i] = m03 + matr[i];\
		m013[i] = m03 + m1[i];\
		m012[i] = m0[i] + m1[i] + matr[i + 4*X/3];\
	}\
	F (t0, vect + 2*X/3, m034);\
	F (t1, vect + X/3, m013);\
	F (t2, vect, m012);\
	F (t3, v1m2, m3);\
	F (t4, v0m2, m0);\
	F (t5, v0m1, m1);\
	for(int i = 0; i < X/3; i++)\
	{\
		rop[i] = t0[i] + t3[i] + t4[i];\
		rop[i + X/3] = t1[i] - t3[i] + t5[i];\
		rop[i + 2*X/3] = t2[i] - t4[i] - t5[i];\
	}\
}

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

void randpoly(poly P)
{
	// Function that generates a random polynomial with all coefficients < RHO.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = (randomint64() % RHO) * (1 + (rand() & 1) * -2);
}


#if !(N % 2)
#define NSBSPLIT N/2
#else
#define NSBSPLIT N/3
#endif
void schoolbookX(__int128* rop, const int64_t* vect, const int64_t* matr)
SCHOOLBOOK(NSBSPLIT)
#undef NSBSPLIT

void toeplitz_vm(__int128* rop, const int64_t* vect, const int64_t* matr)
	#if !(N % 2)
		TOEP22TOP(N, schoolbookX)
	#else
		TOEP33TOP(N, schoolbookX)
	#endif

void pmns_mod_mult_ext_red(__int128* R, const poly A, const poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda  or a binomial as a polynomial used for reduction.
	// Result in R.
	// TODO: add code for other shapes of E.
	
	#if N > 7 && (!(N % 2) || !(N % 3))
		int64_t matr[2*N - 1];
		
		for(int i = 0; i < N-1; i++)
		{
			#ifdef BINOMIAL_A
				matr[i + N - 1] = BINOMIAL_A * B->t[i];
			#else
				matr[i + N - 1] = B->t[i];
			#endif
			matr[i] = B->t[1 + i] * LAMBDA;
		}
		#ifdef BINOMIAL_A
			matr[2*N - 2] = B->t[N - 1] * BINOMIAL_A;
		#else
			matr[2*N - 2] = B->t[N - 1];
		#endif
		toeplitz_vm(R, A->t, matr);
	#else
		#ifdef BINOMIAL_A
			int64_t atab[N];
			for(int i = 0; i < N; i++)
				atab[i] = A->t[i]*BINOMIAL_A;
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


void g_pmns_matvec_mult(__int128* R, const uint64_t* A)
{
	// Function that multiplies A by the short basis G. Result in R.
	
	__int128 somme;
	for(int i = 0; i < N; i++)
	{
		somme = (__int128)G[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (__int128)G[j][i]*A[j];
		
		R[i] += somme;
	}
}

void g1_pmns_matvec_mult(uint64_t* R, const __int128* A)
{
	// Function that multiplies A by G1 = -(G^-1) mod PHI. Result in R.
	
	int64_t somme;
	for(int i = 0; i < N; i++)
	{
		somme = (int64_t)G1[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (int64_t)G1[j][i]*A[j];
		R[i] = somme;
	}
}

void pmns_montg_int_red(poly res, __int128* R)
{
	// Internal reduction of R via the Montgomery-like method.
	uint64_t T[N];
	
	g1_pmns_matvec_mult(T, R);
	
	g_pmns_matvec_mult(R, T);
	
	for(int i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64);
}

void pmns_montg_mult(poly res, const poly A, const poly B)
{
	// Function that multiplies A by B using the Montgomery approach in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	#ifdef BINOMIAL_A
	
	__int128 R[N];
	
	#else
	
	__int128 R[N] = {0};
	
	#endif
	
	pmns_mod_mult_ext_red(R, A, B);
	
	pmns_montg_int_red(res, R);
}

void pmns_rtl_sqandmult(poly res, const poly base, const mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Right to Left version.
	// This algorithm is intended to be constant time.
	
	register uint64_t aux;
	
	poly tmp, num, stok;
	init_polys(N, &tmp, &num, &stok, 0);
	
	poly_copy(res, (const poly) &__theta__);
	poly_copy(stok, base);
	
	for(uint16_t i = 0; i < exponent->deg; i++)
	{
		aux = exponent->t[i];
		for(uint8_t j = 0; j < 32; j++)
		{
			if(aux & 1)
				pmns_montg_mult(tmp, res, stok);
			else
				pmns_montg_mult(tmp, res, (const poly) &__theta__);
			pmns_montg_mult(num, stok, stok);
			aux >>= 1;
			if(aux & 1)
				pmns_montg_mult(res, tmp, num);
			else
				pmns_montg_mult(res, tmp, (const poly) &__theta__);
			pmns_montg_mult(stok, num, num);
			aux >>= 1;
		}
	}
	
	free_polys(tmp, num, stok, 0);
}

void pmns_ltr_sqandmult(poly res, const poly base, const mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Left to Right version.
	// This algorithm is intended to be constant time.
	
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, (const poly) &__theta__);
	
	for(uint16_t i = 0; i < exponent->deg; i++)
	{
		aux = exponent->t[exponent->deg - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			pmns_montg_mult(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				pmns_montg_mult(res, tmp, base);
			else
				pmns_montg_mult(res, tmp, (const poly) &__theta__);
		}
	}
	
	free_poly(tmp);
}

void pmns_montg_ladder(poly res, const poly base, const mpnum exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p.
	// This algorithm is intended to be constant time.
	
	register uint64_t aux, b;
	
	poly tmp, R[2];
	init_poly(N, &tmp);
	R[0] = res;
	R[1] = tmp;
	
	poly_copy(res, (const poly) &__theta__);
	poly_copy(tmp, base);
	
	for(uint16_t i = 0; i < exponent->deg; i++)
	{
		aux = exponent->t[exponent->deg - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			pmns_montg_mult(R[1 - b], R[1 - b], R[b]);
			pmns_montg_mult(R[b], R[b], R[b]);
		}
	}
	
	free_poly(tmp);
}

void convert_binary_to_pmns(poly res, const mpnum op)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (in the Montgomery domain).
	
	uint8_t counter;
	register uint16_t i;
	const uint64_t theta = (1ULL<<THETA);
	__int128 R[N] = {0};
	mpnum stok, tmp;
	init_mpnums(N, &stok, &tmp, 0);
	
	mp_copy(&tmp, op);
	
	mp_mod(&stok, tmp, (const mpnum) &__P__);
	
	for(i = 1; i < N; i++)
		res->t[i] = 0;
	
	// Theta-radix decomposition
	counter = 0;
	res->t[0] = stok->t[0] & (theta - 1);
	for(i = 1; i < stok->deg; i++)
	{
		counter = (counter + THETA) % 64; 
		res->t[i] = (stok->t[i-1] >> counter);
		res->t[i] |= (stok->t[i] << (64 - counter));
		res->t[i] &= (theta - 1);
	}
	counter = (counter + THETA) % 64;
	res->t[i] = (stok->t[i-1] >> counter);
	
	for(i = 0; i < N; i++)
		for(uint16_t j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	pmns_montg_int_red(res, R);
	
	free_mpnums(stok, tmp, 0);
}

void convert_string_to_pmns(poly res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (in the Montgomery domain).
	
	mpnum tmp;
	init_mpnum(N, &tmp);
	
	convert_string_to_binary(&tmp, string);
	
	convert_binary_to_pmns(res, tmp);
	
	free_mpnum(tmp);
}

void convert_pmns_to_binary(mpnum* res, const poly P)
{
	// Function that converts out of the PMNS and into a multiprecision number.
	
	poly a;
	mpnum aux, ag, tmp, tmp2;
	__int128 Quite[N];
	
	init_poly(N, &a);
	init_mpnums(N, &ag, &tmp, &tmp2, 0);
	init_mpnum(1, &aux);
	if((*res)->deg < N)
	{
		free_mpnum(*res);
		init_mpnum(N, res);
	}
	(*res)->deg = N;
	for(uint16_t i = 1; i < N; i++)
	{
		(*res)->t[i] = 0;
		Quite[i] = (__int128) P->t[i];
	}
	Quite[0] = (__int128) P->t[0];
	
	pmns_montg_int_red(a, Quite);
	
	(*res)->sign = 1 - 2 * (a->t[N - 1] < 0);
	(*res)->t[0] = a->t[N - 1] * (*res)->sign;
	for(uint16_t i = 0; i < N - 1; i++)
	{
		mp_mult(&tmp, *res, (const mpnum) &Gi[0]);
		aux->sign = 1 - 2 * (a->t[N - 2 - i] < 0);
		aux->t[0] = a->t[N - 2 - i] * aux->sign;
		mp_add(&tmp2, tmp, aux);
		mp_mod(res, tmp2, (const mpnum) &__P__);
	}
	
	free_mpnums(aux, ag, tmp, tmp2, 0);
	free_poly(a);
}

void pmns_sqandmult(poly res, const char* base, const char* exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Uses Left to Right square and multiply.
	
	poly pbase;
	mpnum mpexponent;
	init_poly(N, &pbase);
	init_mpnum(1, &mpexponent);
	
	convert_string_to_pmns(pbase, base);
	convert_string_to_binary(&mpexponent, exponent);
	
	pmns_ltr_sqandmult(res, pbase, mpexponent);
	
	free_poly(pbase);
	free_mpnum(mpexponent);
}

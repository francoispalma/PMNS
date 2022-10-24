#include <stdio.h>
#include <time.h>

#include "pmns.h"
#include "utilitymp.h"

#define AMNS_MONTG_MULT UNROLLED_amns_montg_mult

inline void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
		
		R[i] = R[i] * LAMBDA;
		
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
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) (A[i + j]) * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) (A[j]) * M[i - j];
	}
}

static inline void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
		for(j = 1; j < N - i; j++)
			R[i] += ((uint64_t)A[i + j]) * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += ((uint64_t)A[j]) * M1[i - j];
	}
}

#endif

#ifdef M_or_B_is_B

static inline void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
			R[i] += (__int128) A[j] * B[j][i];
	}
}

static inline void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Vector-Matrix multiplication between A and B1, result in R.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
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

static inline void UNROLLED_mns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Unrolled version
	uint64_t V[N], V2[N];
	int64_t T[N] = {0};
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	UNROLLED_m1_or_b1_mns_mod_mult_ext_red(T, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = T[i];
	
	UNROLLED_m_or_b_mns_mod_mult_ext_red(R, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = V2[i] + (R[i] >> 64) + (V[i] != 0);
}

inline void UNROLLED_amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Unrolled version
	
	__int128 R[N] = {0};
	
	UNROLLED_mns_mod_mult_ext_red(R, A, B);

	UNROLLED_mns_montg_int_red(res, R);
}

void amns_rtl_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Note that the exponent is a multiprecision
	// number and not a polynomial despite using the poly structure.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		for(j = 0; j < 64; j++)
			{
				if(aux & (1ULL << j))
					AMNS_MONTG_MULT(tmp, res, base);
				else
					AMNS_MONTG_MULT(tmp, res, &__theta__);
				AMNS_MONTG_MULT(res, tmp, tmp);
			}
	}
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			AMNS_MONTG_MULT(tmp, res, base);
		else
			AMNS_MONTG_MULT(tmp, res, &__theta__);
		AMNS_MONTG_MULT(res, tmp, tmp);
		aux >>= 1;
	}
	
	free_poly(tmp);
}

void amns_ltr_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Note that the exponent is a multiprecision
	// number and not a polynomial despite using the poly structure.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux) + 1; j < 64; j++)
	{
		AMNS_MONTG_MULT(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			AMNS_MONTG_MULT(res, tmp, base);
		else
			AMNS_MONTG_MULT(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			AMNS_MONTG_MULT(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				AMNS_MONTG_MULT(res, tmp, base);
			else
				AMNS_MONTG_MULT(res, tmp, &__theta__);
		}
	}
	
	free_poly(tmp);
}

void amns_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Note that the exponent is a multiprecision
	// number and not a polynomial despite using the poly structure.
	
	amns_rtl_sqandmult(res, base, exponent);
}

void amns_montg_ladder(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p. Note that the exponent is a multiprecision
	// number and not a polynomial despite using the poly structure.
	
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
	for(j = __builtin_clz(aux) + 1; j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		AMNS_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
		AMNS_MONTG_MULT(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			AMNS_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
			AMNS_MONTG_MULT(R[b], R[b], R[b]);
		}
	}
	
	free_poly(tmp);
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
	poly stok, tmp;
	init_polys(N, &stok, &tmp, NULL);
	
	if(string[0] == '-')
	{
		printf("ERROR: negative numbers not handled.");
		exit(1);
	}
	
	convert_string_to_multipre(&tmp, string);
	
	mp_utmod(&stok, tmp, &__P__);
	
	if(stok->deg > N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	for(i = 1; i < N; i++)
		res->t[i] = 0;
	
	res->t[0] = ((uint64_t) stok->t[0]) & (rho - 1);
	counter = 0;
	for(i = 1; i < N; i++)
	{
		counter = (counter + 64 - RHO);
		res->t[i] += ((((uint64_t) stok->t[i]) << counter) | (((uint64_t) stok->t[i - 1]) >> (64 - counter))) & (rho - 1);
		if(counter > RHO && i < N - 1)
		{
			counter = counter - RHO;
			res->t[i + 1] = ((((uint64_t) stok->t[i]) << counter) | (((uint64_t) stok->t[i - 1]) >> (64 - counter))) & (rho - 1);
		}
	}
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);

end:
	free_poly(stok);
}

void convert_amns_to_multipre(restrict poly* res, const restrict poly P)
{
	// Function that converts out of the AMNS system and into a multiprecision
	// number. Note that we use the poly structure but res is not a polynomial.
	
	register uint16_t i;
	poly a, aux, ag, tmp;
	__int128 Quite[N];
	
	init_polys(N, &a, &ag, &tmp, NULL);
	init_poly(1, &aux);
	if((*res)->deg < N)
	{
		free_poly(*res);
		init_poly(N, res);
	}
	(*res)->deg = N;
	for(i = 1; i < N; i++)
	{
		(*res)->t[i] = 0;
		Quite[i] = (__int128) P->t[i];
	}
	Quite[0] = (__int128) P->t[0];
	
	mns_montg_int_red(a, Quite);
	
	(*res)->t[0] = a->t[0];
	for(i = 1; i < N; i++)
	{
		aux->t[0] = a->t[i];
		mp_mult(&ag, aux, &Gi[i - 1]);
		
		mp_copy(&tmp, *res);
		mp_add(res, tmp, ag);
	}
	
	mp_copy(&tmp, *res);
	mp_utmod(res, tmp, &__P__);
	
	free_polys(a, aux, ag, tmp, NULL);
}

void __multchecks__(void)
{
	// Used as a debug tool to see if the PMNS correctly gives us the proper
	// results with a few random values.
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	int64_t seed;
	
	srand((unsigned) (time(&seed)));
	
	for(int i = 0; i < 100; i++)
	{
		randpoly(a);
		randpoly(b);
		print(a);
		print(b);
		amns_montg_mult(c, a, b);
		print(c);
	}
	
	free_polys(a, b, c, NULL);
}

void __sqandmultdemo(void)
{
	// TODO: delete it, used to check point to point process.
	const char a[] = "0xffffffffffff", // "7609d69beaadb6de37a7b36cd193b33b120489bd4298534e830eeeaf9a65b15c12268aea1447f610377ea045afc463fb193a531e46cf70052ee6143d782b27aee363d426ad73085f7c24376b676070214cbf1b69f93fd5fdd70b8c77dd2268cbf3f210366b932c7351d9332608fb294ebb44bc7b17bfa3e115dd06c642670d67",
		b[] = "1", // "78a5cc1a942f42e81aa0dc980d3b6a7f987bf9847fb30f8d17c327283745d1365e6bd495a6b4bc2f6a16dca99668ee8591b5b04d12e15d5c4f5f22aafca94fcc3973e7e7714c84b0fb514f862d3444fd36f82f54ccb6c2f2ac3510d3aaea94953533e511076ba29103afb13ba6387001f1055b400b5abfc5b7cd0cdb4b2ebf2a",
		c[] = "4"; // "606ddc8f63ff63abfacb0ced7fcef39d42f0831e9e84459bbffbc1a04b866aaf572571e876a087dc633ab189b40eb861be2f7e194ac5f24ad7886eeb070028bc91a3970ea828cdd5da8eea173d0a38da1bc072837e5835ff2bd266ebcc4788870c7dca82e5b87cd1a844a5120339d2ef1620f1484888538824c5474fe77014e6";
	
	poly A, B, C, aux;
	init_polys(N, &A, &B, &C, &aux, NULL);
	
	convert_string_to_multipre(&B, b);
	
/*	mp_print(B);*/
	
	printf("A = ");
	convert_string_to_amns(A, a);
	
	print(A);
	
	//amns_sqandmult(C, A, B);
	amns_montg_ladder(C, A, B);
	
	printf("\nC = ");
	print(C);
	
	convert_amns_to_multipre(&aux, C);
	
	printf("\n%s\n\n", c);
	mp_print(aux);
	printf("\n");
	
	free_polys(A, B, C, aux, NULL);
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "hpmns.h"
#include "hparams.h"
#include "utilitymp.h"

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


static int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

static void randpoly(poly P)
{
	// Function that generates a random polynomial with all coefficients < 2^RHO.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64()) * (1 + (rand() & 1) * -2);
}

#if N == 9

void schoolbookfortoep3x3(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	for(int i = 0; i < 1; i++)
	{
		rop[i] = 0;
		for(int j = 0; j < 1; j++)
			rop[i] += (__int128) vect[j] * matr[1 - 1 - j + i];
	}
}

void schoolbook3x3(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	for(int i = 0; i < 3; i++)
	{
		rop[i] = 0;
		for(int j = 0; j < 3; j++)
			rop[i] += (__int128) vect[j] * matr[3 - 1 - j + i];
	}
}

void toeplitz_vm3x3(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
TOEP33TOP(3, schoolbookfortoep3x3)

void toeplitz_vm9x9(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
TOEP33TOP(9, schoolbook3x3)

#endif

#ifdef BINOMIAL_A

#if BINOMIAL_B == 1

#if N == 9

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	int64_t matr[2*N - 1];
	
	for(int i = 0; i < N-1; i++)
	{
		matr[i + N - 1] = BINOMIAL_A * B->t[i];
		matr[i] = B->t[1 + i];
	}
	matr[2*N - 2] = B->t[N - 1] * BINOMIAL_A;
	
	toeplitz_vm9x9(R, A->t, matr);
	
}

#else

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	for(int i = 0; i < N; i++)
	{
		somme = (__int128) A->t[0] * B->t[i];
		for(int j = 1; j < i + 1; j++)
			somme += (__int128) A->t[j] * B->t[i - j];
		R[i] = somme * BINOMIAL_A;
	}
	
	for(int i = 0; i < N - 1; i++)
	{
		for(int j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
	}
	
}

#endif

#else

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	for(int i = 0; i < N; i++)
	{
		somme = (__int128) A->t[0] * B->t[i];
		for(int j = 1; j < i + 1; j++)
			somme += (__int128) A->t[j] * B->t[i - j];
		R[i] = somme * BINOMIAL_A;
	}
	
	for(int i = 0; i < N - 1; i++)
	{
		somme = (__int128) A->t[i + 1] * B->t[N - 1];
		for(int j = 2; j < N - i; j++)
			somme += (__int128) A->t[i + j] * B->t[N - j];
		R[i] += somme * BINOMIAL_B;
	}
	
}

#endif

void m_pmns_mod_mult_ext_red(__int128* restrict R,
	const int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	
	R[0] += (__int128)A[0] * -BINOMIAL_TWOTOVERTWO - A[1];
	R[1] += (__int128)A[1] * BINOMIAL_TWOT + A[2];
	for(int i = 2; i < N - 1; i++)
	{
		R[i] += (__int128)A[i] * -BINOMIAL_TWOT + A[i + 1];
	}
	R[N - 1] += (__int128)A[N - 1] * -BINOMIAL_TWOT + A[0];
}

void m1_pmns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	R[0] = (uint64_t)A[N - 2] * -BINOMIAL_TWOT - (uint64_t)A[N - 1];
	R[1] = (uint64_t)A[0] + (uint64_t)A[N - 1] * BINOMIAL_TWOTOVERTWO;
	for(int i = 2; i < N; i++)
	{
		R[i] = (uint64_t)A[i - 2] * -BINOMIAL_TWOT - (uint64_t)A[i - 1];
	}
}


#else

#if N == 9

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	int64_t matr[2*N - 1];
	
	for(int i = 0; i < N-1; i++)
	{
		matr[i + N - 1] = B->t[i];
		matr[i] = B->t[1 + i] * LAMBDA;
	}
	matr[2*N - 2] = B->t[N - 1];
	
	toeplitz_vm9x9(R, A->t, matr);
	
}

#else

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	for(int i = 0; i < N - 1; i++)
	{
		somme = (__int128) A->t[i + 1] * B->t[N - 1];
		for(int j = 2; j < N - i; j++)
			somme += (__int128) A->t[i + j] * B->t[N - j];
		R[i] = somme * LAMBDA;
	}
		
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

#endif

#ifdef TWOTXMONE

void m_pmns_mod_mult_ext_red(__int128* restrict R,
	const uint64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	
	R[0] += (__int128)A[N - 1] * TWOTLAM - A[0];
	for(int i = 1; i < N - 1; i++)
	{
		R[i] += (__int128)A[i - 1] * TWOT - A[i];
	}
	R[N - 1] += (__int128)A[N - 2] * TWOT - A[N - 1] ;
}

void m1_pmns_mod_mult_ext_red(uint64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	R[0] = (uint64_t)A[0] + (uint64_t)A[N - 1] * TWOTLAM;
	for(int i = 1; i < N - 1; i++)
	{
		R[i] = (uint64_t)A[i - 1] * TWOT + (uint64_t)A[i];
	}
	R[N - 1] = (uint64_t)A[N - 2] * TWOT + (uint64_t)A[N - 1];
}

#else

void m_pmns_mod_mult_ext_red(__int128* restrict R,
	const uint64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	
	R[0] += -(__int128)A[0] * GAMMA + (__int128)A[N - 1] * LAMBDA;
	for(int i = 1; i < N - 1; i++)
	{
		R[i] += A[i - 1] - (__int128)A[i] * GAMMA;
	}
	R[N - 1] += -(__int128)A[N - 1] * GAMMA + A[N - 2];
}

#ifdef HOLLOWM1

void m1_pmns_mod_mult_ext_red(uint64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	for(int i = 0; i < N - 2; i++)
	{
		R[i] = -(uint64_t)A[i + 1] - (uint64_t)A[i + 2] * GAMMA;
	}
	R[N - 2] = (uint64_t)A[0] * GAMMALAMM1 - (uint64_t)A[N - 1];
	R[N - 1] = (uint64_t)A[0] * ONELAMM1 + (uint64_t)A[1] * GAMMALAMM1;
}

#else

void m1_pmns_mod_mult_ext_red(uint64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	uint64_t Z = 0;
	
	for(int j = 0; j < N; j++)
		Z += (uint64_t)A[j] * (uint64_t)lastcol[j];
	
	
	R[N - 1] = Z;
	for(int i = 0; i < N - 1; i++)
	{
		Z *= GAMMA;
		Z -= (uint64_t)A[N - 1 - i];
		R[N - 2 - i] = Z;
	}
}

#endif

#endif

#endif

void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Montgomery method.
	uint64_t T[N];
	register uint16_t i;
	
	m1_pmns_mod_mult_ext_red(T, R);
	
	
	/*_poly dummy;
	dummy.t = T;
	dummy.deg = N;
	
	printf("dummy\n");
	poly_print(&dummy);
	printf("\n");*/
	
	m_pmns_mod_mult_ext_red(R, T);
	
	/*printf("[");
	for(int i = 0; i < N; i++)
	{
		printf("0x"); __print128(R[i]); printf(",");
	}
	printf("]\n");*/
	
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64);
}

void pmns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	#ifdef BINOMIAL_A
	
	__int128 R[N];
	
	#else
	
	__int128 R[N] = {0};
	
	#endif
	
	pmns_mod_mult_ext_red(R, A, B);
	
	/*printf("[");
	for(int i = 0; i < N; i++)
	{
		printf("0x"); __print128(R[i]); printf(",");
	}
	printf("]\n");*/
	
	pmns_montg_int_red(res, R);
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
		pmns_montg_mult(c, a, b);
		poly_print(c);
	}
	
	free_polys(a, b, c, NULL);
}

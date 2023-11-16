#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "hpmns.h"
#include "hparams.h"
#include "utilitymp.h"

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

#define TOEP77BOT {\
	__int128 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11,\
	t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23,\
	t24, t25, t26, t27;\
	int64_t v0m6, v0m5, v0m4, v0m3, v0m2, v0m1,\
	v1m6, v1m5, v1m4, v1m3, v1m2,\
	v2m6, v2m5, v2m4, v2m3,\
	v3m6, v3m5, v3m4,\
	v4m6, v4m5,\
	v5m6,\
	m0123456, m1234567, m2345678,\
	m3456789, m456789a, m56789ab,\
	m6789abc;\
	v0m6 = (vect[0] - vect[6]);\
	v1m6 = (vect[1] - vect[6]);\
	v2m6 = (vect[2] - vect[6]);\
	v3m6 = (vect[3] - vect[6]);\
	v4m6 = (vect[4] - vect[6]);\
	v5m6 = (vect[5] - vect[6]);\
	v0m5 = (vect[0] - vect[5]);\
	v1m5 = (vect[1] - vect[5]);\
	v2m5 = (vect[2] - vect[5]);\
	v3m5 = (vect[3] - vect[5]);\
	v4m5 = (vect[4] - vect[5]);\
	v0m4 = (vect[0] - vect[4]);\
	v1m4 = (vect[1] - vect[4]);\
	v2m4 = (vect[2] - vect[4]);\
	v3m4 = (vect[3] - vect[4]);\
	v0m3 = (vect[0] - vect[3]);\
	v1m3 = (vect[1] - vect[3]);\
	v2m3 = (vect[2] - vect[3]);\
	v0m2 = (vect[0] - vect[2]);\
	v1m2 = (vect[1] - vect[2]);\
	v0m1 = (vect[0] - vect[1]);\
	int64_t m67 = matr[6] + matr[7], m678 = m67 + matr[8], m6789 = m678 + matr[9],\
	m6789a = m6789 + matr[10], m6789ab = m6789a + matr[11];\
	int64_t m45 = matr[4] + matr[5], m345 = matr[3] + m45, m2345 = matr[2] + m345,\
	m12345 = matr[1] + m2345;\
	m0123456 = matr[0] + m12345 + matr[6];\
	m1234567 = m12345 + m67;\
	m2345678 = m2345 + m678;\
	m3456789 = m345 + m6789;\
	m456789a = m45 + m6789a;\
	m56789ab = matr[5] + m6789ab;\
	m6789abc = matr[12] + m6789ab;\
	t0 = (__int128)vect[6] * m0123456;\
	t1 = (__int128)vect[5] * m1234567;\
	t2 = (__int128)vect[4] * m2345678;\
	t3 = (__int128)vect[3] * m3456789;\
	t4 = (__int128)vect[2] * m456789a;\
	t5 = (__int128)vect[1] * m56789ab;\
	t6 = (__int128)vect[0] * m6789abc;\
	t7 = (__int128)v0m6 * matr[6];\
	t8 = (__int128)v1m6 * matr[5];\
	t9 = (__int128)v2m6 * matr[4];\
	t10 = (__int128)v3m6 * matr[3];\
	t11 = (__int128)v4m6 * matr[2];\
	t12 = (__int128)v5m6 * matr[1];\
	t13 = (__int128)v0m5 * matr[7];\
	t14 = (__int128)v1m5 * matr[6];\
	t15 = (__int128)v2m5 * matr[5];\
	t16 = (__int128)v3m5 * matr[4];\
	t17 = (__int128)v4m5 * matr[3];\
	t18 = (__int128)v0m4 * matr[8];\
	t19 = (__int128)v1m4 * matr[7];\
	t20 = (__int128)v2m4 * matr[6];\
	t21 = (__int128)v3m4 * matr[5];\
	t22 = (__int128)v0m3 * matr[9];\
	t23 = (__int128)v1m3 * matr[8];\
	t24 = (__int128)v2m3 * matr[7];\
	t25 = (__int128)v0m2 * matr[10];\
	t26 = (__int128)v1m2 * matr[9];\
	t27 = (__int128)v0m1 * matr[11];\
	rop[0] = t0 + t7 + t8 + t9 + t10 + t11 + t12;\
	rop[1] = t1 - t12 + t13 + t14 + t15 + t16 + t17;\
	rop[2] = t2 - t11 - t17 + t18 + t19 + t20 + t21;\
	rop[3] = t3 - t10 - t16 - t21 + t22 + t23 + t24;\
	rop[4] = t4 - t9 - t15 - t20 - t24 + t25 + t26;\
	rop[5] = t5 - t8 - t14 - t19 - t23 - t26 + t27;\
	rop[6] = t6 - t7 - t13 - t18 - t22 - t25 - t27;\
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

#if N % 3
#define NSBSPLIT N/2
#else
#define NSBSPLIT N/3
#endif
void schoolbookX(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
SCHOOLBOOK(NSBSPLIT)
#undef NSBSPLIT

void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	#if N == 7
		TOEP77BOT
	#elif N % 3
		TOEP22TOP(N, schoolbookX)
	#else
		TOEP33TOP(N, schoolbookX)
	#endif

void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	#if !(N % 2) || !(N % 3)
		int64_t matr[2*N - 1];
		
		for(int i = 0; i < N-1; i++)
		{
		#ifdef BINOMIAL_A
			matr[i + N - 1] = BINOMIAL_A * B->t[i];
		#else
			matr[i + N - 1] = B->t[i];
		#endif
		#if LAMBDA == 1
			matr[i] = B->t[1 + i];
		#else
			matr[i] = B->t[1 + i] * LAMBDA;
		#endif
		}
		#ifdef BINOMIAL_A
			matr[2*N - 2] = B->t[N - 1] * BINOMIAL_A;
		#else
			matr[2*N - 2] = B->t[N - 1];
		#endif
		toeplitz_vm(R, A->t, matr);
	#else
		__int128 somme;
		
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
			for(int i = 0; i < N; i++)
				btab[i] = B->t[i]*LAMBDA;
		#endif//TOEP77BOT
		for(int i = 0; i < N; i++)
		{
			somme = (__int128) atab[0] * B->t[i];
			for(int j = 1; j < i + 1; j++)
				somme += (__int128) atab[j] * B->t[i - j];
			R[i] = somme;
		}
		
		for(int i = 0; i < N - 1; i++)
		{
			for(int j = 1; j < N - i; j++)
				R[i] += (__int128) A->t[i + j] * btab[N - j];
		}
	#endif
}


void m_pmns_mod_mult_ext_red(__int128* restrict R,
	const uint64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	#ifdef BINOMIAL_TWOTOVERTWO
	R[0] += (__int128)A[0] * -BINOMIAL_TWOTOVERTWO - A[1];
	R[1] += (__int128)A[1] * BINOMIAL_TWOT + A[2];
	for(int i = 2; i < N - 1; i++)
	{
		R[i] += (__int128)A[i] * -BINOMIAL_TWOT + A[i + 1];
	}
	R[N - 1] += (__int128)A[N - 1] * -BINOMIAL_TWOT + A[0];
	#else
	#ifdef TWOTXMONE
	R[0] += (__int128)A[N - 1] * TWOTLAM - A[0];
	for(int i = 1; i < N; i++)
	{
		R[i] += (__int128)A[i - 1] * TWOT - A[i];
	}
	#else
	#ifdef BETH
	R[0] += -(__int128)A[0] * GAMMA + (__int128)A[N - 1] * BETH;
	#else
	R[0] += -(__int128)A[0] * GAMMA + (__int128)A[N - 1] * LAMBDA;
	#endif
	for(int i = 1; i < N - 1; i++)
	{
		R[i] += A[i - 1] - (__int128)A[i] * GAMMA;
	}
	#ifdef GIMEL
	R[N - 1] += -(__int128)A[N - 1] * GIMEL + A[N - 2];
	#else
	R[N - 1] += -(__int128)A[N - 1] * GAMMA + A[N - 2];
	#endif
	#endif
	#endif
}

void m1_pmns_mod_mult_ext_red(uint64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	#ifdef BINOMIAL_TWOTOVERTWO
		R[0] = (uint64_t)A[N - 2] * -BINOMIAL_TWOT - (uint64_t)A[N - 1];
		R[1] = (uint64_t)A[0] + (uint64_t)A[N - 1] * BINOMIAL_TWOTOVERTWO;
		for(int i = 2; i < N; i++)
		{
			R[i] = (uint64_t)A[i - 2] * -BINOMIAL_TWOT - (uint64_t)A[i - 1];
		}
	#else
		#ifdef TWOTXMONE
			R[0] = (uint64_t)A[0] + (uint64_t)A[N - 1] * TWOTLAM;
			for(int i = 1; i < N; i++)
			{
				R[i] = (uint64_t)A[i - 1] * TWOT + (uint64_t)A[i];
			}
		#else
			#ifdef HOLLOWM1
				#ifndef GIMEL
				for(int i = 0; i < N - 2; i++)
				{
					R[i] = -(uint64_t)A[i + 1] - (uint64_t)A[i + 2] * GAMMA;
				}
				R[N - 2] = -(uint64_t)A[0] * GAMMALAMM1 - (uint64_t)A[N - 1];
				R[N - 1] = -(uint64_t)A[0] * ONELAMM1 - (uint64_t)A[1] * GAMMALAMM1;
				#else
				for(int i = 0; i < N - 2; i++)
				{
					R[i] = -(uint64_t)A[i + 1] - (uint64_t)A[i + 2] * GAMMA;
				}
				R[N - 2] = -(uint64_t)A[0] * GIMEL - (uint64_t)A[N - 1];
				R[N - 1] = -(uint64_t)A[0] - (uint64_t)A[1] * GAMMA;
				#endif
			#else
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
			#endif
		#endif
	#endif
}

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

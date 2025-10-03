#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <omp.h>

#define mpn_mullo_n __MPN(mullo_n)
__GMP_DECLSPEC void mpn_mullo_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

#define mpn_redc_n __MPN(redc_n)
__GMP_DECLSPEC void mpn_redc_n (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);

#define NTEST 501
#define NSAMPLES 1001
#define NBTHREADZ 6

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))
#define HOW(X) (X>>64)

#include "toeplitzmacros.h"

#include "redparams.h"

__int128 M2[2*N-1];
__int128 BigM[2*N-1];

typedef struct
{
	uint16_t deg;
	int64_t *t;
} _poly, *poly;

typedef struct
{
	uint16_t deg;
	__int128 *t;
} _poly128, *poly128;

static inline void poly128_print(const restrict poly128 P)
{
	int64_t hi; uint64_t lo;
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
	{
		hi = (int64_t)(P->t[i]>>64);
		lo = (uint64_t)(P->t[i]);
		if(hi >= 0)
			printf("0x%lx%016lx, ", hi, lo);
		else
			printf("-0x%lx%016lx, ", -hi - 1, -lo);
	}
	hi = (int64_t)(P->t[P->deg - 1]>>64);
	lo = (uint64_t)(P->t[P->deg - 1]);
	if(hi >= 0)
		printf("0x%lx%016lx]\n", hi, lo);
	else
		printf("-0x%lx%016lx]\n", -hi - 1, -lo);
}

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline __int128 randomint128(void)
{
	return (((__int128)randomint64())<<64) | randomint64();
}

void randpoly(poly P)
{
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->t[i] = randomint64() % RHO;
		if(randomint64() & 1)
			P->t[i] *= -1;
	}
}

void randpoly128(poly128 P)
{
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->t[i] = (((__int128)(randomint64() % RHOhi))<<64) | randomint64();
		if(randomint64() & 1)
			P->t[i] *= -1;
	}
}

void gmp_mult_wrapper(restrict poly res, const restrict poly A, const restrict poly B)
{
	mp_limb_t blank[2*LENGTH_OF_P];
	mp_limb_t tmp[2*LENGTH_OF_P] = {0};
	mp_limb_t mid[16][LENGTH_OF_P] = {0};
	uint64_t* At = (uint64_t*) A->t;
	uint64_t* Bt = (uint64_t*) B->t;
	#pragma omp parallel for num_threads(8)
	for(int i = 0; i < 16; i++)
	{
		if(i == 0)
			mpn_mul_n(mid[i], At, Bt, LENGTH_OF_P/4);
		else if(i == 1)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt, LENGTH_OF_P/4);
		else if(i == 2)
			mpn_mul_n(mid[i], At, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 3)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt, LENGTH_OF_P/4);
		else if(i == 4)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 5)
			mpn_mul_n(mid[i], At, Bt + LENGTH_OF_P/2, LENGTH_OF_P/4);
		else if(i == 6)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt, LENGTH_OF_P/4);
		else if(i == 7)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 8)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt + LENGTH_OF_P/2, LENGTH_OF_P/4);
		else if(i == 9)
			mpn_mul_n(mid[i], At, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 10)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 11)
			mpn_mul_n(mid[i], At + 2*LENGTH_OF_P/4, Bt + 2*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 12)
			mpn_mul_n(mid[i], At + 1*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 13)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + 2*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 14)
			mpn_mul_n(mid[i], At + 2*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 15)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
	}
	int chac = 1;
	mpn_add_n(tmp, tmp, mid[0], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[1], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[2], LENGTH_OF_P/2);
	chac = 2;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[3], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[4], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[5], LENGTH_OF_P/2);
	chac = 3;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[6], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[7], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[8], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[9], LENGTH_OF_P/2);
	chac = 4;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[10], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[11], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[12], LENGTH_OF_P/2);
	chac = 5;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[13], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[14], LENGTH_OF_P/2);
	mpn_add_n(tmp+3*LENGTH_OF_P/2, tmp+3*LENGTH_OF_P/2, mid[15], LENGTH_OF_P/2);
	mpn_tdiv_qr(blank, (uint64_t*)res->t, 0, tmp, 2*LENGTH_OF_P, __P__, LENGTH_OF_P);
}

void gmp_montg_wrapper(restrict poly res, const restrict poly A, const restrict poly B)
{
	mp_limb_t tmp[2*LENGTH_OF_P];
	mp_limb_t mid[16][LENGTH_OF_P] = {0};
	uint64_t* At = (uint64_t*) A->t;
	uint64_t* Bt = (uint64_t*) B->t;
	#pragma omp parallel for num_threads(8)
	for(int i = 0; i < 16; i++)
	{
		if(i == 0)
			mpn_mul_n(mid[i], At, Bt, LENGTH_OF_P/4);
		else if(i == 1)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt, LENGTH_OF_P/4);
		else if(i == 2)
			mpn_mul_n(mid[i], At, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 3)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt, LENGTH_OF_P/4);
		else if(i == 4)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 5)
			mpn_mul_n(mid[i], At, Bt + LENGTH_OF_P/2, LENGTH_OF_P/4);
		else if(i == 6)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt, LENGTH_OF_P/4);
		else if(i == 7)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 8)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/4, Bt + LENGTH_OF_P/2, LENGTH_OF_P/4);
		else if(i == 9)
			mpn_mul_n(mid[i], At, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 10)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 11)
			mpn_mul_n(mid[i], At + 2*LENGTH_OF_P/4, Bt + 2*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 12)
			mpn_mul_n(mid[i], At + 1*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 13)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + 2*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 14)
			mpn_mul_n(mid[i], At + 2*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
		else if(i == 15)
			mpn_mul_n(mid[i], At + 3*LENGTH_OF_P/4, Bt + 3*LENGTH_OF_P/4, LENGTH_OF_P/4);
	}
	int chac = 1;
	mpn_add_n(tmp, tmp, mid[0], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[1], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[2], LENGTH_OF_P/2);
	chac = 2;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[3], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[4], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[5], LENGTH_OF_P/2);
	chac = 3;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[6], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[7], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[8], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[9], LENGTH_OF_P/2);
	chac = 4;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[10], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[11], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[12], LENGTH_OF_P/2);
	chac = 5;
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[13], LENGTH_OF_P/2);
	tmp[(chac+2)*LENGTH_OF_P/4] += mpn_add_n(tmp + chac*LENGTH_OF_P/4, tmp + chac*LENGTH_OF_P/4, mid[14], LENGTH_OF_P/2);
	mpn_add_n(tmp+3*LENGTH_OF_P/2, tmp+3*LENGTH_OF_P/2, mid[15], LENGTH_OF_P/2);
	/*#pragma omp parallel for num_threads(2)
	for(int i = 0; i < 4; i++)
	{
		if(i == 0)
			mpn_mul_n(mid[i], At, Bt, LENGTH_OF_P/2);
		else if(i == 1)
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt, LENGTH_OF_P/2);
		else if(i == 2)
			mpn_mul_n(mid[i], At, Bt + LENGTH_OF_P/2, LENGTH_OF_P/2);
		else
			mpn_mul_n(mid[i], At + LENGTH_OF_P/2, Bt + LENGTH_OF_P/2, LENGTH_OF_P/2);
	}
	tmp[3*LENGTH_OF_P/2] = mpn_add_n(tmp + LENGTH_OF_P/2, mid[1], mid[2], LENGTH_OF_P);
	tmp[LENGTH_OF_P] += mpn_add_n(tmp, tmp, mid[0], LENGTH_OF_P);
	mpn_add_n(tmp + LENGTH_OF_P, tmp + LENGTH_OF_P, mid[3], LENGTH_OF_P);*/
	//mpn_mul_n(tmp, At, Bt, LENGTH_OF_P);
	mpn_redc_n((uint64_t*)res->t, tmp, __P__, LENGTH_OF_P, __invP__);
}

static inline void pmns_mod_mult_ext_red(__int128* restrict res, const restrict poly A, const restrict poly B)
{
	__int128 Res[N];
	int64_t matr[2*N-1];
	
	// We construct the Toeplitz Matrices
	for(int i = 0; i < N; i++)
		matr[i + N - 1] = B->t[i];
	for(int i = 0; i < N - 1; i++)
		matr[i] = B->t[1 + i] * LAMBDA;
	
	// Res <- A * B mod E
	toeplitz_vm(Res, A->t, matr);
	
	for(int i = 0; i < N; i++)
		res[i] = Res[i];
}


void pmns_montg_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t T[N];
	__int128 Res[N];
	__int128 aux[N];
	
	pmns_mod_mult_ext_red(Res, A, B);
	
	/*printf("Res = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(Res[i]>>64), (uint64_t)Res[i]);
	printf("]\n");*/
	
	// T <- Res * M1 mod E mod PHI
	m1toeplitz_vm(T, Res);
	
	/*printf("T = [");
	for(int i = 0; i < N; i++)
		printf("%ld, ", T[i]);
	printf("]\n");*/
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	/*for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			Res[i] += (__int128) T[j] * M[N - 1 - j + i];
	}*/
	
		
	/*printf("aux = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(Res[i]>>64), (uint64_t)Res[i]);
	printf("]\n");*/
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		res->t[i] = ((__int128)Res[i] + aux[i]) >> 64;
}

/*void toeplitz_vm3x3_leveldown(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/3)

void toeplitz_vm3x3_multithread(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[N/3], t1[N/3], t2[N/3], t3[N/3], t4[N/3], t5[N/3];
	
	int64_t m03, v1m2[N/3], v0m2[N/3], v0m1[N/3], m034[N-1], m013[N-1], m012[N-1], m0[N-1], m1[N-1], m3[N-1];
	
	
	for(int i = 0; i < N/3; i++)
	{
		v1m2[i] = vect[N/3 + i] - vect[2*N/3 + i];
		v0m1[i] = vect[i] - vect[N/3 + i];
		v0m2[i] = vect[i] - vect[2*N/3 + i];
	}
	
	for(int i = 0; i < N-1; i++)
	{
		m0[i] = matr[i + 2*N/3];
		m1[i] = matr[i + N];
		m3[i] = matr[i + N/3];
		m03 = m0[i] + m3[i];
		m034[i] = m03 + matr[i];
		m013[i] = m03 + m1[i];
		m012[i] = m0[i] + m1[i] + matr[i + 4*N/3];
	}
	
	#pragma omp parallel for num_threads(6)
	for(int i = 0; i < 6; i++)
	{
		if(i == 0)
			mtoep10(t0, vect + 2*N/3, m034);
		else if(i == 1)
			mtoep10(t1, vect + N/3, m013);
		else if(i == 2)
			mtoep10(t2, vect, m012);
		else if(i == 3)
			mtoep10(t3, v1m2, m3);
		else if(i == 4)
			mtoep10(t4, v0m2, m0);
		else if(i == 5)
			mtoep10(t5, v0m1, m1);
	}
	/*#pragma omp parallel sections num_threads(6)
	{
		#pragma omp section
		mtoep10(t0, vect + 2*N/3, m034);
		#pragma omp section
		mtoep10(t1, vect + N/3, m013);
		#pragma omp section
		mtoep10(t2, vect, m012);
		#pragma omp section
		mtoep10(t3, v1m2, m3);
		#pragma omp section
		mtoep10(t4, v0m2, m0);
		#pragma omp section
		mtoep10(t5, v0m1, m1);
	}*//*
	
	for(int i = 0; i < N/3; i++)
	{
		rop[i] = t0[i] + t3[i] + t4[i];
		rop[i + N/3] = t1[i] - t3[i] + t5[i];
		rop[i + 2*N/3] = t2[i] - t4[i] - t5[i];
	}
}

void mtoeplitz_vm3x3_multithread(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[N/3], t1[N/3], t2[N/3], t3[N/3], t4[N/3], t5[N/3],v1m2[N/3], v0m2[N/3], v0m1[N/3];
	
	int64_t m03, m034[N-1], m013[N-1], m012[N-1], m0[N-1], m1[N-1], m3[N-1];
	
	
	for(int i = 0; i < N/3; i++)
	{
		v1m2[i] = (__int128)vect[N/3 + i] - vect[2*N/3 + i];
		v0m1[i] = (__int128)vect[i] - vect[N/3 + i];
		v0m2[i] = (__int128)vect[i] - vect[2*N/3 + i];
	}
	
	for(int i = 0; i < N-1; i++)
	{
		m0[i] = matr[i + 2*N/3];
		m1[i] = matr[i + N];
		m3[i] = matr[i + N/3];
		m03 = m0[i] + m3[i];
		m034[i] = m03 + matr[i];
		m013[i] = m03 + m1[i];
		m012[i] = m0[i] + m1[i] + matr[i + 4*N/3];
	}
	
	#pragma omp parallel for num_threads(6)
	for(int i = 0; i < 6; i++)
	{
		if(i == 0)
			mmtoep10(t0, vect + 2*N/3, m034);
		else if(i == 1)
			mmtoep10(t1, vect + N/3, m013);
		else if(i == 2)
			mmtoep10(t2, vect, m012);
		else if(i == 3)
			mmmtoep10(t3, v1m2, m3);
		else if(i == 4)
			mmmtoep10(t4, v0m2, m0);
		else if(i == 5)
			mmmtoep10(t5, v0m1, m1);
	}
	
	for(int i = 0; i < N/3; i++)
	{
		rop[i] = t0[i] + t3[i] + t4[i];
		rop[i + N/3] = t1[i] - t3[i] + t5[i];
		rop[i + 2*N/3] = t2[i] - t4[i] - t5[i];
	}
}

void pmns_montg_mult_multithread(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t T[N];
	__int128 Res[N];
	
	{
		int64_t matrix[2*N - 1];
		
		for(int i = 0; i < N-1; i++)
		{
			matrix[i + N - 1] = B->t[i];
			matrix[i] = B->t[1 + i] * LAMBDA;
		}
		matrix[2*N - 2] = B->t[N - 1];
		
		/*#pragma omp parallel for shared(A, matrix, Res) private(somme) num_threads(NBTHREADZ)
		for(int i = 0; i < N; i++)
		{
			somme = 0;
			for(int j = 0; j < N; j++)
				somme += (__int128) A->t[j] * matrix[N - 1 - j + i];
			Res[i] = somme;
		}*//*
		toeplitz_vm3x3_multithread(Res, A->t, matrix);
	}
	
	// T <- Res * M1 mod E mod PHI
	{
		uint64_t somme;
		#pragma omp parallel for shared(T, M1, Res) private(somme) num_threads(NBTHREADZ)
		for(int i = 0; i < N; i++)
		{
			somme = (uint64_t)M1[N-1+i]*Res[0];
			for (int j = 1; j < N; j++)
				somme += (uint64_t)M1[N-1-j+i]*Res[j];
			T[i] = somme;
		}
	}
	__int128 aux[N];
	// aux <- T * M
	{
		/*__int128 somme;
		#pragma omp parallel for shared(T, M, Res) private(somme) num_threads(NBTHREADZ)
		for(int i = 0; i < N; i++)
		{
			somme = (__int128)M[N-1+i]*T[0];
			for (int j = 1; j < N; j++)
				somme += (__int128)M[N-1-j+i]*T[j];
			
			Res[i] += somme;
		}*//*
		
		mtoeplitz_vm3x3_multithread(aux, T, M);
	}
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		res->t[i] = ((__int128)Res[i]+aux[i]) >> 64;
}*/

void pmns_montg_mult128(restrict poly128 res, const restrict poly128 A, const restrict poly128 B)
{
	__int128 T[N];
	int64_t Thi[N];
	uint64_t Tlo[N];
	__int128 Reshihi[N];
	__int128 Reslohi[N];
	__int128 Reshilo[N] = {0};
	__int128 Reslolo[N];
	__int128 Reslo[N];
	__int128 auxhihi[N];
	__int128 auxlohi[N];
	__int128 auxhilo[N] = {0};
	
	int64_t Ahi[N];
	uint64_t Alo[N];
	for(int i = 0; i < N; i++)
	{
		Ahi[i] = (int64_t)(A->t[i]>>64);
		Alo[i] = (uint64_t)(A->t[i]);
	}
	
	int64_t mathi[2*N-1];
	uint64_t matlo[2*N-1];
	__int128 bigmat;
	for(int i = 0; i < N; i++)
	{
		mathi[N-1+i] = (int64_t)(B->t[i]>>64);
		matlo[N-1+i] = (uint64_t)(B->t[i]);
	}
	for(int i = 0; i < N-1; i++)
	{
		bigmat = (__int128)B->t[i+1]*LAMBDA;
		mathi[i] = (int64_t)(bigmat>>64);
		matlo[i] = (uint64_t)(bigmat);
	}
	
	toeplitz_vm(Reshihi, Ahi, mathi);
	/*toeplitz_vm(Reshilo, Ahi, matlo);
	mtoeptop(Reslohi, Alo, mathi);
	mtoeptop(Reslolo, Alo, matlo);*/
	for(int i = 0; i < N; i++)
	{
		Reslohi[i] = 0;
		for(int j = 0; j < N; j++)
		{
			Reslohi[i] += (__int128)Alo[j] * mathi[N-1-j+i] + (__int128)Ahi[j] * matlo[N-1-j+i];
		}
	}
	for(int i = 0; i < N; i++)
	{
		Reslolo[i] = 0;
		for(int j = 0; j < N; j++)
		{
			Reslolo[i] += (__int128)Alo[j] * matlo[N-1-j+i];
		}
	}
	
	/*_poly128 fordisplay;
	__int128 displaytab[N];
	fordisplay.deg = N;/*
	fordisplay.t = Reslolo;
	printf("Reslolo = ZZX("); poly128_print(&fordisplay); printf(")\n");
	fordisplay.t = Reslohi;
	printf("Reslohi = ZZX("); poly128_print(&fordisplay); printf(")\n");
	fordisplay.t = Reshihi;
	printf("Reshihi = ZZX("); poly128_print(&fordisplay); printf(")\n");/**/
	
	for(int i = 0; i < N; i++)
		Reslo[i] = Reslolo[i] + (((__int128)Reshilo[i] + Reslohi[i])<<64);
	/*fordisplay.t = Reslolo;
	printf("Reslo = ZZX("); poly128_print(&fordisplay); printf(")\n");/**/
	
	// T <- Res * M1 mod E mod PHI
	toep128(T, Reslo, M2);
	/*fordisplay.t = M2 + N-1;
	printf("M2 = ZZX("); poly128_print(&fordisplay); printf(")\n");
	fordisplay.t = T;
	printf("T = ZZX("); poly128_print(&fordisplay); printf(")\n");/**/
	
	for(int i = 0; i < N; i++)
	{
		Thi[i] = (int64_t)(T[i]>>64);
		Tlo[i] = (uint64_t)(T[i]);
	}
	
	__int128 auxlolo[N];
	// aux <- T * M
	toeplitz_vm(auxhihi, Thi, Mhi);
	/*toeplitz_vm(auxhilo, Thi, Mlo);
	mtoeptop(auxlohi, Tlo, Mhi);/**/
	//toeplitz_ll(auxlolo, Tlo, Mlo);
	
	for(int i = 0; i < N; i++)
	{
		auxlohi[i] = 0;
		for(int j = 0; j < N; j++)
		{
			auxlohi[i] += (__int128)Tlo[j] * Mhi[N-1-j+i] + (__int128)Thi[j] * Mlo[N-1-j+i];
		}
	}
	
	for(int i = 0; i < N; i++)
	{
		auxlolo[i] = 0;
		for(int j = 0; j < N; j++)
		{
			auxlolo[i] += (__int128)Tlo[j] * Mlo[N-1-j+i];
		}
	}
	
	/*fordisplay.t = auxhihi;
	printf("auxhihi = ZZX("); poly128_print(&fordisplay); printf(")\n");
	fordisplay.t = auxlohi;
	printf("auxlohi = ZZX("); poly128_print(&fordisplay); printf(")\n");/**/
	//fordisplay.t = auxlolo;
	//printf("auxlolo = ZZX("); poly128_print(&fordisplay); printf(")\n");
	
	// res <- (Res + aux)/PHI
	__int128 tmp[N];
	for(int i = 0; i < N; i++)
	{
		tmp[i] = ((__int128)Reshilo[i] + Reslohi[i] + auxhilo[i] + auxlohi[i]);
		res->t[i] = ((__int128)Reshihi[i] + auxhihi[i]) + (tmp[i]>>64);
		tmp[i] = ((__int128)auxhilo[i] + auxlohi[i]);
	}
	//fordisplay.t = tmp;
	//printf("tmp = ZZX("); poly128_print(&fordisplay); printf(")\n");
}

static inline _Bool add_overflow(unsigned __int128* restrict a, const unsigned __int128 b)
{
	const unsigned __int128 tmp = *a;
	*a += b;
	return *a < tmp;
}

void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	int64_t tmplo;
	uint64_t Shi[N];
	unsigned __int128 tmp, aux, Slo[N], Tlo[N] = {0}, Thi[N] = {0};
	
	for(int i = 0; i < N - 1; i++)
	{
		Slo[i] = (unsigned __int128) LOW(A->t[i + 1]) * LOW(B->t[N - 1]);
		Shi[i] = 0;
		for(int j = 2; j < N - i; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) LOW(A->t[i + j]) * LOW(B->t[N - j]));
		}
		
		aux = (unsigned __int128) LOW(Slo[i]) * LAMBDA;
		tmp = (unsigned __int128) HI(Slo[i]) * LAMBDA + HI(aux);
		Slo[i] = ((unsigned __int128) tmp << 64) | LOW(aux);
		Shi[i] = (unsigned __int128) Shi[i] * LAMBDA + HI(tmp);
		
		for(int j = 0; j < i + 1; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) LOW(A->t[j]) * LOW(B->t[i - j]));
		}
	}
	Slo[N - 1] = (unsigned __int128) LOW(A->t[0]) * LOW(B->t[N - 1]);
	Shi[N - 1] = add_overflow(Slo + N - 1, (unsigned __int128) LOW(A->t[1]) * LOW(B->t[N - 2]));
	for(int j = 2; j < N; j++)
	{
		Shi[N - 1] += add_overflow(Slo + N - 1, (unsigned __int128) LOW(A->t[j]) * LOW(B->t[N - 1 - j]));
	}
	
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * M2lo[N-1+i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * M2lo[N-1-j+i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) Mlo[N-1-i+j] * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) Mhi[N-1-i+j] * tmplo + HIGH(tmp);
		}
	}
	
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HI(Slo[i]) + ((__int128)Shi[i] << 64) + HIGH(Tlo[i]);
		Slo[i] = 0;
		for(int j = 1; j < N - i; j++)
		{
			Slo[i] += (__int128) LOW(A->t[i + j]) * HIGH(B->t[N - j]) +
				(__int128) HIGH(A->t[i + j]) * LOW(B->t[N - j]);
		}
		
		Slo[i] *= LAMBDA;
		
		for(int j = 0; j < i + 1; j++)
		{
			Slo[i] += (__int128) HIGH(A->t[j]) * LOW(B->t[i - j]) +
				(__int128) LOW(A->t[j]) * HIGH(B->t[i - j]);
		}
		Slo[i] += Thi[i] + 1;
		Thi[i] = 0;
		Tlo[i] = 0;
	}
	
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * M2lo[N-1+i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * M2lo[N-1-j+i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) Mlo[N-1-i+j] * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) Mhi[N-1-i+j] * tmplo + HIGH(tmp);
		}
	}
	
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HIGH(Slo[i]) + HIGH(Tlo[i]) + 1;
		tmp = 0;
		for(int j = 1; j < N - i; j++)
		{
			tmp += (__int128) HIGH(A->t[i + j]) * HIGH(B->t[N - j]);
		}
		
		tmp *= LAMBDA;
		
		for(int j = 0; j < i + 1; j++)
		{
			tmp += (__int128) HIGH(A->t[j]) * HIGH(B->t[i - j]);
		}
		Thi[i] += tmp;
	}
	
	for(int i = 0; i < N; i++)
	{
		res->t[i] = Thi[i];
	}
}

void pmns_barrett_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t T[N];
	__int128 Res[N];
	int64_t divi[N];
	__int128 m1aux[N];
	__int128 aux[N];
	
	pmns_mod_mult_ext_red(Res, A, B);
	
	/*printf("Res = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(Res[i]>>64), (uint64_t)Res[i]);
	printf("]\n");/**/
	
	for(int i = 0; i < N; i++)
	{
		divi[i] = (int64_t)((Res[i]+(1ULL<<DIVBYRHO))>>(DIVBYRHO+1));
	}
	/*printf("divi = [");
	for(int i = 0; i < N; i++)
		printf("%ld, ", divi[i]);
	printf("]\n");/**/
	
	toeplitz_vm(m1aux, divi, Gprime);
	/*printf("m1aux = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(m1aux[i]>>64), (uint64_t)m1aux[i]);
	printf("]\n");/**/
	for(int i = 0; i < N; i++)
	{
		/*tmp = (int64_t)(m1aux[i]>>63);
		T[i] = (tmp>>1) + (tmp&1);/**/
		T[i] = (int64_t)((m1aux[i] + (1ULL<<63))>>64);
	}
	/*printf("T = [");
	for(int i = 0; i < N; i++)
		printf("%ld, ", T[i]);
	printf("]\n");/**/
	// aux <- T * M
	toeplitz_vm(aux, T, M);
	/*printf("aux = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(aux[i]>>64), (uint64_t)aux[i]);
	printf("]\n");/**/
	
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		res->t[i] = (int64_t)((__int128)Res[i] - aux[i]);
}

void pmns_plant_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	__int128 Res[N];
	__int128 aux[N];
	int64_t hi[N];
	/*int64_t tmphi;
	__int128 tmpmid;*/
	
	pmns_mod_mult_ext_red(Res, A, B);
	
	/*printf("Res = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(Res[i]>>64), (uint64_t)Res[i]);
	printf("]\n");*/
	
	/*for(int i = 0; i < N; i++)
	{
		Reshi[i] = (int64_t)(Res[i]>>64);
		Reslo[i] = ((int64_t)Res[i]);
	}*/
	
	// T <- -Res * M1 mod E mod PHI^2
	/*m1toep20(hi, Reshi, M2hi);
	toeplitz_vm(mid1, Reshi, M2lo);
	toeplitz_vm(mid2, Reslo, M2hi);
	toeplitz_vm(lo, Reslo, M2lo);*/
	
	/*for(int i = 0; i < N; i++)
	{
		mid1[i] += (int64_t)(lo[i]>>64);
		hi[i] += (int64_t)((mid1[i] + mid2[i])>>64) + 1;
	}*/
		
	/*for(int i = 0; i < N; i++)
	{
		//tmp[i] = 0;
		hi[i] = 1;
		for(int j = 0; j < N; j++)
			hi[i] += (int64_t)((((__int128) Res[j] * M2[N-1-j+i]))>>64);
	}*/
			//tmp[i] += (__int128) Res[j] * M2[N-1-j+i];
		/*{
			tmphi = (Res[j]>>64) * M2hi[N-1-j+i];
			tmpmid = (__int128)M2hi[N-1-j+i] * ((int64_t)Res[j]) + (__int128)M2lo[N-1-j+i] * (Res[j]>>64);
			tmp[i] += tmpmid + ((__int128)tmphi<<64);
		}*/
			//
		//hi[i] = (int64_t)((tmp[i]+(1ULL<<63))>>64) + 1;
	
	__int128 tmp[N];
	toep128(tmp, Res, M2);
	for(int i = 0; i < N; i++)
		hi[i] = (int64_t)((tmp[i])>>64) + 1;
	
	/*printf("T = [");
	for(int i = 0; i < N; i++)
		printf("%ld, ", T[i]);
	printf("]\n");*/
	
	// aux <- (T/PHI + (1,...,1)) * M mod E
	//toeplitz_vm(aux, hi, M);
	
	for(int i = 0; i< N; i++)
	{
		aux[i] = 0;
		for(int j = 0; j < N; j++)
			aux[i] += (__int128)hi[j] * M[N-1-j+i];
	}
	
	/*for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			Res[i] += (__int128) T[j] * M[N - 1 - j + i];
	}*/
	
		
	/*printf("aux = [");
	for(int i = 0; i < N; i++)
		printf("0x%lx%016lx, ", (int64_t)(Res[i]>>64), (uint64_t)Res[i]);
	printf("]\n");*/
	
	// res <- (aux)/PHI
	for(int i = 0; i < N; i++)
		res->t[i] = ((__int128)aux[i]+(1ULL<<63)) >> 64;
}

static inline void poly_print(const restrict poly P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[P->deg - 1]);
}



/**** Measurements procedures according to INTEL white paper

	"How to benchmark code execution times on INTEL IA-32 and IA-64"

*****/

void quicksort(uint64_t* t, int n)
{
	if (n > 0)
	{
	/* partitionning */
	int i, j, temp;
	
	j=0;
	for (i=0; i<n-1; i++)
	{
		/* at least as big as the pivot */
		if (t[i] < t[n-1])
		{
		temp = t[j];
		t[j] = t[i];
		t[i] = temp;
		j++;
		}
	}
	/*replacing the pivot */
	temp = t[j];
	t[j] = t[n-1];
	t[n-1] = temp;
	
	quicksort(t, j);
	quicksort(&t[j+1], n-j-1);
	}
}

uint64_t *quartiles(uint64_t *tab, uint64_t size)
{
	uint64_t *result ;
	uint64_t aux ;
	
	result = malloc(3*sizeof(uint64_t));
	quicksort(tab, size);
	aux = size >> 2;
	if (size % 4) aux++;
	// Q1
	result[0] = tab[aux-1];
	// Mediane
	// size is odd hence it's easy
	result[1]	= tab[(size+1)/2 - 1];
	// Q3
	aux = (3*size) >> 2;
	if ((3*size) % 4) aux++;
	result[2]	= tab[aux - 1];
	
	return result;
}

static inline uint64_t cpucyclesStart(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n    "
			"RDTSC\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

static inline uint64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			"CPUID\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_poly a, b, c;
	int64_t atab[N], btab[N], ctab[N];
	a.deg = N; b.deg = N; c.deg = N;
	a.t = atab;
	b.t = btab;
	c.t = ctab;
	poly A = &a, B = &b, C = &c;
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(A);
		randpoly(B);
		pmns_mult(C, A, B);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(A);
		randpoly(B);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			printf("\b%d\t%d\r", i, j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				pmns_mult(C, A, B);
				pmns_mult(B, C, A);
				pmns_mult(A, B, C);
			}
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			if(timermin > diff_t) timermin = diff_t;
			else if(timermax < diff_t) timermax = diff_t;
			cycles[j]=diff_t;
		}
		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles,NTEST);
		medianTimer += statTimer[1];
		free(statTimer);
	}
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_bench128(void (*pmns_mult)(restrict poly128, const restrict poly128, const restrict poly128), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_poly128 a, b, c;
	__int128 atab[N], btab[N], ctab[N];
	a.deg = N; b.deg = N; c.deg = N;
	a.t = atab;
	b.t = btab;
	c.t = ctab;
	poly128 A = &a, B = &b, C = &c;
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly128(A);
		randpoly128(B);
		pmns_mult(C, A, B);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly128(A);
		randpoly128(B);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			printf("\b%d\t%d\r", i, j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				pmns_mult(C, A, B);
				pmns_mult(B, C, A);
				pmns_mult(A, B, C);
			}
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			if(timermin > diff_t) timermin = diff_t;
			else if(timermax < diff_t) timermax = diff_t;
			cycles[j]=diff_t;
		}
		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles,NTEST);
		medianTimer += statTimer[1];
		free(statTimer);
	}
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_expbench(void (*pmns_exp)(restrict poly, const restrict poly, const int64_t* restrict), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc((NTEST/10 + 1),sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_poly a, b, c;
	int64_t atab[N], btab[N], ctab[N];
	a.deg = N; b.deg = N; c.deg = N;
	a.t = atab;
	b.t = btab;
	c.t = ctab;
	poly A = &a, B = &b, C = &c;
	int64_t exponent[LENGTH_OF_P];
	int indexfortestlater;
	
	for(int i=0;i<(NTEST/10 + 1);i++)
	{
	// Here we "heat" the cache memory.
		randpoly(A);
		for(int j = 0; j< LENGTH_OF_P; j++)
			exponent[j] = randomint64();
		indexfortestlater = LENGTH_OF_P - 1;
		while(indexfortestlater >= 0 && exponent[indexfortestlater] == __P__[indexfortestlater])
			indexfortestlater--;
		if(indexfortestlater >= 0 && exponent[indexfortestlater] > __P__[indexfortestlater])
			for(int j = 0; j < LENGTH_OF_P; j++)
				exponent[j] -= __P__[j];
		
		pmns_exp(B, A, exponent);
	}
	
	for(int i=0;i<NSAMPLES/10;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(A);
		for(int j = 0; j< LENGTH_OF_P; j++)
			exponent[j] = randomint64();
		indexfortestlater = LENGTH_OF_P - 1;
		while(indexfortestlater >= 0 && exponent[indexfortestlater] == __P__[indexfortestlater])
			indexfortestlater--;
		if(indexfortestlater >= 0 && exponent[indexfortestlater] > __P__[indexfortestlater])
			for(int j = 0; j < LENGTH_OF_P; j++)
				exponent[j] -= __P__[j];
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,(NTEST/10 + 1)*sizeof(uint64_t));
		for(int j=0;j<(NTEST/10 + 1);j++)
		{
			printf("\b%d\t%d\r", i, j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				pmns_exp(B, A, exponent);
				pmns_exp(C, B, exponent);
				pmns_exp(A, C, exponent);
			}
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			if(timermin > diff_t) timermin = diff_t;
			else if(timermax < diff_t) timermax = diff_t;
			cycles[j]=diff_t;
		}
		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles,(NTEST/10 + 1));
		medianTimer += statTimer[1];
		free(statTimer);
	}
	
	printf("                                          \r");
	
	free(cycles);
	return meanTimermin/(NSAMPLES/10)/W; // We divide by W since we measured W calls.
}

/****

	End of section.

*****/


/*void gmp_ltr_wrapper(poly res, const restrict poly base, const int64_t* restrict exponent)
{
	uint64_t aux;
	__int128 tmp[LENGTH_OF_P*2];
	__int128 blank[LENGTH_OF_P*2];
	for(int i = 0; i < LENGTH_OF_P; i++)
		res->t[i] = base->t[i];
	
	for(uint16_t i = 0; i < LENGTH_OF_P; i++)
	{
		aux = exponent[LENGTH_OF_P - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			mpn_mul_n(tmp, res->t, res->t, LENGTH_OF_P);
			mpn_tdiv_qr(blank, res->t, 0, tmp, LENGTH_OF_P*2, __P__, LENGTH_OF_P);
			if(aux & (1ULL << (63 - j)))
			{
				mpn_mul_n(tmp, res->t, base->t, LENGTH_OF_P);
				mpn_tdiv_qr(blank, res->t, 0, tmp, LENGTH_OF_P*2, __P__, LENGTH_OF_P);
			}
		}
	}
}*/

void pmns_ltr_montg(poly res, const restrict poly base, const int64_t* restrict exponent)
{
	uint64_t aux;
	
	for(int i = 0; i < N; i++)
		res->t[i] = base->t[i];
	
	for(uint16_t i = 0; i < LENGTH_OF_P; i++)
	{
		aux = exponent[LENGTH_OF_P - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			pmns_montg_mult(res, res, res);
			if(aux & (1ULL << (63 - j)))
				pmns_montg_mult(res, res, base);
		}
	}
}

void pmns_ltr_barrett(poly res, const restrict poly base, const int64_t* restrict exponent)
{
	uint64_t aux;
	
	for(int i = 0; i < N; i++)
		res->t[i] = base->t[i];
	
	for(uint16_t i = 0; i < LENGTH_OF_P; i++)
	{
		aux = exponent[LENGTH_OF_P - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			pmns_barrett_mult(res, res, res);
			if(aux & (1ULL << (63 - j)))
				pmns_barrett_mult(res, res, base);
		}
	}
}

void pmns_ltr_plant(poly res, const restrict poly base, const int64_t* restrict exponent)
{
	uint64_t aux;
	
	for(int i = 0; i < N; i++)
		res->t[i] = base->t[i];
	
	__int128 bigmat[2*N-1];
	__int128 bigres[N];
	int64_t hi[N];
	__int128 second[N];
	
	bigmattoep(bigmat+N-1, base->t, M2);
	for(int i = 0; i < N-1; i++)
		bigmat[i] = (__int128)bigmat[N+i]*LAMBDA;
	
	
	for(uint16_t i = 0; i < LENGTH_OF_P; i++)
	{
		aux = exponent[LENGTH_OF_P - 1 - i];
		for(uint8_t j = 0; j < 64; j++)
		{
			pmns_plant_mult(res, res, res);
			if(aux & (1ULL << (63 - j)))
			{
				bigmattoep(bigres, res->t, bigmat);
				for(int k = 0; k < N; k++)
					hi[i] = (int64_t)((bigres[k])>>64) + 1;
				for(int k = 0; k < N; k++)
				{
					second[k] = 0;
					for(int l = 0; l < N; l++)
						second[l] += (__int128)hi[l] * M[N-1-l+k];
				}
			}
		}
	}
}



int main(void)
{
	for(int i = 0; i < 2*N - 1; i++)
		M2[i] = (((__int128)M2hi[i]<<64) | (M2lo[i]));
	for(int i = 0; i < 2*N - 1; i++)
		BigM[i] = (((__int128)Mhi[i]<<64) | (Mlo[i]));
	
	time_t seed;
	srand((unsigned) (time(&seed)));
	
#if N > 48
	int W = 3;
#else
	int W = 33;
#endif
	
	uint64_t pcycles, mcycles = 0, bcycles = 0;
#ifndef NOBENCH
	//printf("cycles %ld\n", do_bench(gmp_mult_wrapper,W));
	//printf("cycles %ld\n", do_bench(gmp_montg_wrapper,W));
	//printf("Plantard-like cycles %ld\n", do_bench(pmns_plant_mult, W));
#ifndef EXPBENCH
	pcycles = do_bench(pmns_plant_mult, W);
#endif
#ifndef NOMONTG
	//printf("Montgomery-like cycles %ld\n", do_bench(pmns_montg_mult,W));
	//printf("cycles %ld\n", do_bench(pmns_montg_mult_multithread,W));
	mcycles = do_bench(pmns_montg_mult, W);
#ifndef NOBARRETT
	//printf("Barrett-like cycles %ld\n", do_bench(pmns_barrett_mult,W));
	bcycles = do_bench(pmns_barrett_mult, W);
#endif
#endif
	//printf("cycles128 %ld\n", do_bench128(amns128_montg_mult,W));
	
	//printf("gnump cycles %ld\n", do_expbench(gmp_ltr_wrapper,3));
#ifdef EXPBENCH
	//printf("plantard cycles %ld\n", do_expbench(pmns_ltr_plant,3));
	//printf("montgom cycles %ld\n", do_expbench(pmns_ltr_montg,3));
	//printf("barrett cycles %ld\n", do_expbench(pmns_ltr_barrett,3));
	pcycles = do_expbench(pmns_ltr_plant,3)/1000;
	mcycles = do_expbench(pmns_ltr_montg,3)/1000;
	bcycles = do_expbench(pmns_ltr_barrett,3)/1000;
#endif
	char mstr[10], bstr[10];
	sprintf(mstr, mcycles ? "%ld" : "-", mcycles);
	sprintf(bstr, bcycles ? "%ld" : "-", bcycles);
	printf("||\t%d\t||\t%ld\t||\t%s\t||\t%s\t||\n", N, pcycles, mstr, bstr); 
#endif
	/*_poly A, B, C;
	int64_t atab[N], btab[N], ctab[N];
	A.deg = N; B.deg = N; C.deg = N;
	A.t = atab;
	B.t = btab;
	C.t = ctab;
	randpoly(&A);randpoly(&B);
	pmns_plant_mult(&C,&A,&B);
	printf("A = ZZX("); poly_print(&A);printf(")\n");
	printf("B = ZZX("); poly_print(&B);printf(")\n");
	printf("C = ZZX("); poly_print(&C);printf(")\n");/**/
	
	/*_poly128 A, B, C;
	__int128 atab[N], btab[N], ctab[N];
	A.deg = N; B.deg = N; C.deg = N;
	A.t = atab;
	B.t = btab;
	C.t = ctab;
	randpoly128(&A);randpoly128(&B);
	pmns_montg_mult128(&C,&A,&B);
	printf("A = ZZX("); poly128_print(&A);printf(")\n");
	printf("B = ZZX("); poly128_print(&B);printf(")\n");
	printf("C = ZZX("); poly128_print(&C);printf(")\n");
	amns128_montg_mult(&C, &A, &B);
	printf("A = ZZX("); poly128_print(&A);printf(")\n");
	printf("B = ZZX("); poly128_print(&B);printf(")\n");
	printf("C = ZZX("); poly128_print(&C);printf(")\n");*/
	
	/**/FILE *fpointer = fopen("mlog", "w+");
	fclose(fpointer);
	fpointer = freopen("mlog", "a+", stdout);
	
	_poly A, B, C;
	int64_t atab[N], btab[N], ctab[N];
	A.deg = N; B.deg = N; C.deg = N;
	A.t = atab;
	B.t = btab;
	C.t = ctab;
	for(int i = 0; i < 1000; i++)
	{
		randpoly(&A);
		randpoly(&B);
		poly_print(&A);
		poly_print(&B);
		pmns_montg_mult(&C, &A, &B);
		//pmns_montg_mult_multithread(&C, &A, &B);
		poly_print(&C);
	}
	
	fclose(fpointer);
	fpointer = fopen("blog", "w+");
	fclose(fpointer);
	fpointer = freopen("blog", "a+", stdout);
	
	for(int i = 0; i < 1000; i++)
	{
		randpoly(&A);
		randpoly(&B);
		poly_print(&A);
		poly_print(&B);
		pmns_barrett_mult(&C, &A, &B);
		poly_print(&C);
	}
	
	fclose(fpointer);
	fpointer = fopen("plog", "w+");
	fclose(fpointer);
	fpointer = freopen("plog", "a+", stdout);
	
	for(int i = 0; i < 1000; i++)
	{
		randpoly(&A);
		randpoly(&B);
		poly_print(&A);
		poly_print(&B);
		pmns_plant_mult(&C, &A, &B);
		poly_print(&C);
	}
	
	fclose(fpointer);/**/
	
	/*FILE *fpointer = fopen("plog", "w+");
	fclose(fpointer);
	fpointer = freopen("plog", "a+", stdout);
	
	_poly128 A, B, C;
	__int128 atab[N], btab[N], ctab[N];
	A.deg = N; B.deg = N; C.deg = N;
	A.t = atab;
	B.t = btab;
	C.t = ctab;
	for(int i = 0; i < 1000; i++)
	{
		randpoly128(&A);
		randpoly128(&B);
		poly128_print(&A);
		poly128_print(&B);
		amns128_montg_mult(&C, &A, &B);
		poly128_print(&C);
	}
	
	fclose(fpointer);*/
	
	return 0;
}

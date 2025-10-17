#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include <string.h>

#define NTEST 501
#define NSAMPLES 1001

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))
#define HOW(X) (X>>64)

#include "toeplitzmacros.h"
#include "equparams.h"

__int128 M2[2*N-1];
__int128 trans_T[N];

typedef struct
{
	uint16_t deg;
	int64_t *t;
} _poly, *poly;

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
	
	// T <- Res * M1 mod E mod PHI
	m1toeplitz_vm(T, Res);
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		res->t[i] = ((__int128)Res[i] + aux[i]) >> 64;
}

void pmns_plant_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	__int128 Res[N];
	__int128 aux[N];
	int64_t hi[N];
	
	pmns_mod_mult_ext_red(Res, A, B);
	
	__int128 tmp[N];
	toep128(tmp, Res, M2);
	for(int i = 0; i < N; i++)
		hi[i] = (int64_t)((tmp[i])>>64) + 1;
	
	
	for(int i = 0; i< N; i++)
	{
		aux[i] = 0;
		for(int j = 0; j < N; j++)
			aux[i] += (__int128)hi[j] * M[N-1-j+i];
	}
	
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

uint64_t do_equalitybench(_Bool (*pmns_equtest)(const restrict poly, const restrict poly), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_poly a, b;
	uint8_t cccchek = 0;
	int64_t atab[N], btab[N];
	a.deg = N; b.deg = N;
	a.t = atab;
	b.t = btab;
	poly A = &a, B = &b;
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(A);
		randpoly(B);
		cccchek = cccchek + pmns_equtest(A, B);
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
			for(uint64_t soak=0; soak < W/2; soak++)
			{
				cccchek = cccchek + pmns_equtest(A, B) + pmns_equtest(B, A);
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
	
	printf("%d", cccchek);
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_doubleequalitybench(_Bool (*pmns_equtest)(const __int128* restrict, const __int128* restrict), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	uint8_t cccchek = 0;
	__int128 A[N], B[N];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		for(int k = 0; k < N; k++)
		{
			A[k] = randomint128() % WRHOCARRE;
			B[k] = randomint128() % WRHOCARRE;
		}
		cccchek = cccchek + pmns_equtest(A, B);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		for(int k = 0; k < N; k++)
		{
			A[k] = randomint128();
			B[k] = randomint128();
		}
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			printf("\b%d\t%d\r", i, j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/2; soak++)
			{
				cccchek = cccchek + pmns_equtest(A, B) + pmns_equtest(B, A);
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
	
	printf("%d", cccchek);
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

/****

	End of section.

*****/

_Bool naive_equality_test(const restrict poly A, const restrict poly B)
{
	int64_t ctab[N]; _poly C; C.deg = N; C.t = ctab;
	// C  = A - B
	for(int i = 0; i < N; i++)
		C.t[i] = A->t[i] - B->t[i];
	
	// Conversion of C to binary
	uint64_t c_check[LENGTH_OF_P + 1] = {0};
	uint64_t tmp[LENGTH_OF_P + 1] = {0};
	uint64_t soak[N+1], res[LENGTH_OF_P+1];
	c_check[0] = C.t[0];
	for(int i = 1; i < N; i++)
	{
		tmp[LENGTH_OF_P] = mpn_mul_1(tmp, __Gi__[i], LENGTH_OF_P, C.t[i]);
		mpn_add_n(c_check, c_check, tmp, LENGTH_OF_P + 1);
	}
	mpn_tdiv_qr(soak, res, 0, c_check, LENGTH_OF_P + 1, __P__, LENGTH_OF_P);
	
	// Check that the result is 0
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < LENGTH_OF_P))
	{
		flag = res[i] == 0;
		i++;
	}
	return flag;
}

_Bool naive_equality_doubletest(const __int128 * restrict A, const __int128 * restrict B)
{
	int64_t ctab[N]; _poly C; C.deg = N; C.t = ctab;
	int64_t T[N];
	__int128 Res[N];
	__int128 aux[N];
	
	for(int i = 0; i < N; i++)
		Res[i] = A[i] - B[i];
	// T <- Res * M1 mod E mod PHI
	m1toeplitz_vm(T, Res);
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	// C <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		C.t[i] = ((__int128)Res[i] + aux[i]) >> 64;
	
	// Conversion of C to binary
	uint64_t c_check[LENGTH_OF_P + 1] = {0};
	uint64_t tmp[LENGTH_OF_P + 1] = {0};
	uint64_t soak[N+1], res[LENGTH_OF_P+1];
	c_check[0] = C.t[0];
	for(int i = 1; i < N; i++)
	{
		tmp[LENGTH_OF_P] = mpn_mul_1(tmp, __Gi__[i], LENGTH_OF_P, C.t[i]);
		mpn_add_n(c_check, c_check, tmp, LENGTH_OF_P + 1);
	}
	mpn_tdiv_qr(soak, res, 0, c_check, LENGTH_OF_P + 1, __P__, LENGTH_OF_P);
	
	// Check that the result is 0
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < LENGTH_OF_P))
	{
		flag = res[i] == 0;
		i++;
	}
	return flag;
}

_Bool barrett_equality_test(const restrict poly A, const restrict poly B)
{
	int64_t T[N];
	int64_t S[N];
	__int128 tmp[N];
	int64_t ctab[N]; _poly C; C.deg = N; C.t = ctab;
	// C  = A - B
	for(int i = 0; i < N; i++)
		C.t[i] = A->t[i] - B->t[i];
	
	// T <- round((C * G')/PHI)
	toeplitz_vm(tmp, C.t, Gprime);
	for(int i = 0; i < N; i++)
		T[i] = (tmp[i] + (1ULL<<63))>>64;
	
	// S <- T * G
	m1toep20(S, T, M);
	
	// Check that result is correct
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = S[i] == C.t[i];
		i++;
	}
	return flag;
}

_Bool barrett_equality_doubletest(const __int128 * restrict A, const __int128 * restrict B)
{
	int64_t T[N];
	int64_t S[N];
	__int128 tmp[N];
	int64_t ctab[N]; _poly C; C.deg = N; C.t = ctab;
	__int128 Res[N];
	__int128 aux[N];
	
	for(int i = 0; i < N; i++)
		Res[i] = A[i] - B[i];
	// T <- Res * M1 mod E mod PHI
	m1toeplitz_vm(T, Res);
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	// C <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		C.t[i] = ((__int128)Res[i] + aux[i]) >> 64;
	
	// T <- round((C * G')/PHI)
	toeplitz_vm(tmp, C.t, Gprime);
	for(int i = 0; i < N; i++)
		T[i] = (tmp[i] + (1ULL<<63))>>64;
	
	// S <- T * G
	m1toep20(S, T, M);
	
	// Check that result is correct
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = S[i] == C.t[i];
		i++;
	}
	return flag;
}

#ifdef TRANSLATIONPOSSIBLE
_Bool translation_equality_halftest(const restrict poly A, const restrict poly B)
{
	int64_t check[N];
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = trans_T[i] + A->t[i] - B->t[i];
	
	int64_t T[N];
	__int128 aux[N];
	
	// T <- C * M1 mod E mod PHI
	m1toeplitz_vm(T, C);
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		check[i] = ((__int128)C[i] + aux[i]) >> 64;
	
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = check[i] == 0;
		i++;
	}
	return flag;
}

_Bool translation_equality_test(const __int128 * restrict A, const __int128 * restrict B)
{
	int64_t check[N];
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = trans_T[i] + A[i] - B[i];
	
	int64_t T[N];
	__int128 aux[N];
	
	// T <- C * M1 mod E mod PHI
	m1toeplitz_vm(T, C);
	
	// aux <- T * M
	mtoeplitz_vm(aux, T);
	
	// res <- (Res + aux)/PHI
	for(int i = 0; i < N; i++)
		check[i] = ((__int128)C[i] + aux[i]) >> 64;
	
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = check[i] == 0;
		i++;
	}
	return flag;
}
#endif

_Bool plantard_equality_halftest(const restrict poly A, const restrict poly B)
{
	int64_t S[N];
	int64_t C[N];
	__int128 aux[N];
	int64_t hi[N];
	__int128 tmp[N];
	// C  = A - B
	for(int i = 0; i < N; i++)
		C[i] = A->t[i] - B->t[i];
	
	// hi <- (C*M1 mod E mod PHI^2)/PHI + 1
	bigmattoep(tmp, C, M2);
	for(int i = 0; i < N; i++)
		hi[i] = (int64_t)((tmp[i])>>64) + 1;
	
	// aux <- hi * M mod E
	for(int i = 0; i< N; i++)
	{
		aux[i] = 0;
		for(int j = 0; j < N; j++)
			aux[i] += (__int128)hi[j] * M[N-1-j+i];
	}
	
	// S <- round(aux/PHI)
	for(int i = 0; i < N; i++)
		S[i] = ((__int128)aux[i]+(1ULL<<63)) >> 64;
	
	// Check that result is correct
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = S[i] == 0;
		i++;
	}
	return flag;
}

_Bool plantard_equality_test(const __int128 * restrict A, const __int128 * restrict B)
{
	int64_t S[N];
	__int128 C[N];
	__int128 aux[N];
	int64_t hi[N];
	__int128 tmp[N];
	// C  = A - B
	for(int i = 0; i < N; i++)
		C[i] = A[i] - B[i];
	
	// hi <- (C*M1 mod E mod PHI^2)/PHI + 1
	toep128(tmp, C, M2);
	for(int i = 0; i < N; i++)
		hi[i] = (int64_t)((tmp[i])>>64) + 1;
	
	// aux <- hi * M mod E
	for(int i = 0; i< N; i++)
	{
		aux[i] = 0;
		for(int j = 0; j < N; j++)
			aux[i] += (__int128)hi[j] * M[N-1-j+i];
	}
	
	// S <- round(aux/PHI)
	for(int i = 0; i < N; i++)
		S[i] = ((__int128)aux[i]+(1ULL<<63)) >> 64;
	
	// Check that result is correct
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		flag = S[i] == 0;
		i++;
	}
	return flag;
}

int main(void)
{
	for(int i = 0; i < 2*N - 1; i++)
		M2[i] = (((__int128)M2hi[i]<<64) | (M2lo[i]));
	#ifdef TRANSLATIONPOSSIBLE
	for(int i = 0; i < N; i++)
		trans_T[i] = (((__int128)trans_Thi[i]<<64) | (trans_Tlo[i]));
	#endif
	
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	_poly A, B;
	int64_t atab[N], btab[N];
	A.deg = N; B.deg = N;
	A.t = atab; B.t = btab;
	
	randpoly(&A); randpoly(&B);
	
	printf("naive\t%s", naive_equality_test(&A, &B) ? "True" : "False");
	printf("\t%s", naive_equality_test(&A, &A) ? "True" : "False");
	printf("\t%s\n", naive_equality_test(&B, &B) ? "True" : "False");
	printf("barrett\t%s", barrett_equality_test(&A, &B) ? "True" : "False");
	printf("\t%s", barrett_equality_test(&A, &A) ? "True" : "False");
	printf("\t%s\n", barrett_equality_test(&B, &B) ? "True" : "False");
	#ifdef TRANSLATIONPOSSIBLE
	printf("transla\t%s", translation_equality_halftest(&A, &B) ? "True" : "False");
	printf("\t%s", translation_equality_halftest(&A, &A) ? "True" : "False");
	printf("\t%s\n", translation_equality_halftest(&B, &B) ? "True" : "False");
	#endif
	printf("plantar\t%s", plantard_equality_halftest(&A, &B) ? "True" : "False");
	printf("\t%s", plantard_equality_halftest(&A, &A) ? "True" : "False");
	printf("\t%s\n", plantard_equality_halftest(&B, &B) ? "True" : "False");
	
#if N > 47
	int W = 6;
#else
	int W = 30;
#endif
	
	printf("Montgomery-like cycles %ld\n", do_bench(pmns_montg_mult,W));
	printf("naive test cycles %ld\n", do_equalitybench(naive_equality_test,W));
	printf("barrett test cycles %ld\n", do_equalitybench(barrett_equality_test,W));
	#ifdef TRANSLATIONPOSSIBLE
	printf("translat test cycles %ld\n", do_equalitybench(translation_equality_halftest,W));
	#endif
	printf("plantard test cycles %ld\n", do_equalitybench(plantard_equality_halftest,W));
	
	__int128 dA[N], dB[N];
	for(int i = 0; i < N; i++)
	{
		dA[i] = randomint128() % WRHOCARRE;
		dB[i] = randomint128() % WRHOCARRE;
	}
	
	printf("double naive\t%s", naive_equality_doubletest(dA, dB) ? "True" : "False");
	printf("\t%s", naive_equality_doubletest(dA, dA) ? "True" : "False");
	printf("\t%s\n", naive_equality_doubletest(dB, dB) ? "True" : "False");
	printf("double barrett\t%s", barrett_equality_doubletest(dA, dB) ? "True" : "False");
	printf("\t%s", barrett_equality_doubletest(dA, dA) ? "True" : "False");
	printf("\t%s\n", barrett_equality_doubletest(dB, dB) ? "True" : "False");
	#ifdef TRANSLATIONPOSSIBLE
	printf("double transla\t%s", translation_equality_test(dA, dB) ? "True" : "False");
	printf("\t%s", translation_equality_test(dA, dA) ? "True" : "False");
	printf("\t%s\n", translation_equality_test(dB, dB) ? "True" : "False");
	#endif
	printf("double plantar\t%s", plantard_equality_test(dA, dB) ? "True" : "False");
	printf("\t%s", plantard_equality_test(dA, dA) ? "True" : "False");
	printf("\t%s\n", plantard_equality_test(dB, dB) ? "True" : "False");
	
	printf("naive doubletest cycles %ld\n", do_doubleequalitybench(naive_equality_doubletest,W));
	printf("barrett doubletest cycles %ld\n", do_doubleequalitybench(barrett_equality_doubletest,W));
	#ifdef TRANSLATIONPOSSIBLE
	printf("translat doubletest cycles %ld\n", do_doubleequalitybench(translation_equality_test,W));
	#endif
	printf("plantard doubletest cycles %ld\n", do_doubleequalitybench(plantard_equality_test,W));
	
	return 0;
}


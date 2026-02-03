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

#include "toeplitzmacros.h"
#include "equparams.h"


__int128 M2[2*N-1];
__int128 G2[N][N];
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
	
	
	mtoeplitz_vm(aux, hi);
	
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
	uint64_t counter = 0;
	uint64_t cccchek = 0;
	int64_t atab[N], btab[N];
	a.deg = N; b.deg = N;
	a.t = atab;
	b.t = btab;
	poly A = &a, B = &b;
	#ifdef EQUALITY
	int8_t dummy;
	#endif
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(A);
		randpoly(B);
		#ifdef EQUALITY
		for(int j = 0; j < N; j++)
			B->t[j] = A->t[j];
		for(int i2 = 0; i2 < N; i2++)
		{
			dummy = (rand()%7) - 3;
			for(int j2 = 0; j2 < N; j2++)
				B->t[j2] += dummy*G[i2][j2];
		}
		#endif
		cccchek = cccchek + pmns_equtest(A, B);
		#ifdef EQUALITY
		counter++;
		#endif
		
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(A);
		randpoly(B);
		#ifdef EQUALITY
		for(int j = 0; j < N; j++)
			B->t[j] = A->t[j];
		for(int i2 = 0; i2 < N; i2++)
		{
			dummy = (rand()%7) - 3;
			for(int j2 = 0; j2 < N; j2++)
				B->t[j2] += dummy*G[i2][j2];
		}
		#endif
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
			#ifdef EQUALITY
			for(uint64_t soak=0; soak < W/2; soak++)
			{
				counter++;
				counter++;
			}
			#endif
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
	if(cccchek != counter)
		printf("\nError:%ld\t%ld\n", cccchek, counter);
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_doubleequalitybench(_Bool (*pmns_equtest)(const __int128* restrict, const __int128* restrict), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	uint64_t counter = 0;
	uint64_t cccchek = 0;
	__int128 A[N], B[N];
	#ifdef EQUALITY
	int64_t dummy;
	#endif
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		for(int k = 0; k < N; k++)
			B[k] = 0;
		for(int k = 0; k < N; k++)
		{
			A[k] = randomint128() % WRHOCARRE;
			#ifdef EQUALITY
			B[k] += A[k];
			dummy = rand()-(1ULL<<31);
			for(int j2 = 0; j2 < N; j2++)
				B[j2] += (__int128)dummy*G[k][j2];
			#else
			B[k] = randomint128() % WRHOCARRE;
			#endif
		}
		cccchek = cccchek + pmns_equtest(A, B);
		#ifdef EQUALITY
		counter++;
		#endif
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		for(int k = 0; k < N; k++)
			B[k] = 0;
		for(int k = 0; k < N; k++)
		{
			A[k] = randomint128() % WRHOCARRE;
			#ifdef EQUALITY
			B[k] += A[k];
			dummy = rand()-(1ULL<<31);
			for(int j2 = 0; j2 < N; j2++)
				B[j2] += (__int128)dummy*G[k][j2];
			#else
			B[k] = randomint128() % WRHOCARRE;
			#endif
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
			#ifdef EQUALITY
			for(uint64_t soak=0; soak < W/2; soak++)
			{
				counter++;
				counter++;
			}
			#endif
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
	if(cccchek != counter)
		printf("\nError:%ld\t%ld\n", cccchek, counter);
	
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
	uint64_t cpos[LENGTH_OF_P + 1] = {0};
	uint64_t cneg[LENGTH_OF_P + 1] = {0};
	uint64_t tmp[LENGTH_OF_P + 1] = {0};
	uint64_t soak[N+1], res[LENGTH_OF_P+1];
	if(C.t[0] > 0)
		cpos[0] = C.t[0];
	else
		cneg[0] = -C.t[0];
	for(int i = 1; i < N; i++)
	{
		tmp[LENGTH_OF_P] = mpn_mul_1(tmp, __Gi__[i-1], LENGTH_OF_P, C.t[i]*(C.t[i]>0) - C.t[i]*(C.t[i]<0));
		if(C.t[i] > 0)
		{
			mpn_add_n(cpos, cpos, tmp, LENGTH_OF_P + 1);
		}
		else
		{
			mpn_add_n(cneg, cneg, tmp, LENGTH_OF_P + 1);
		}
	}
	if(mpn_cmp(cpos, cneg, LENGTH_OF_P + 1)>=0)
		mpn_sub_n(c_check, cpos, cneg, LENGTH_OF_P + 1);
	else
		mpn_sub_n(c_check, cneg, cpos, LENGTH_OF_P + 1);
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
	
	uint64_t c_check[LENGTH_OF_P + 1] = {0};
	uint64_t cpos[LENGTH_OF_P + 1] = {0};
	uint64_t cneg[LENGTH_OF_P + 1] = {0};
	uint64_t tmp[LENGTH_OF_P + 1] = {0};
	uint64_t soak[N+1], res[LENGTH_OF_P+1];
	if(C.t[0] > 0)
		cpos[0] = C.t[0];
	else
		cneg[0] = -C.t[0];
	for(int i = 1; i < N; i++)
	{
		tmp[LENGTH_OF_P] = mpn_mul_1(tmp, __Gi__[i-1], LENGTH_OF_P, C.t[i]*(C.t[i]>0) - C.t[i]*(C.t[i]<0));
		if(C.t[i] > 0)
		{
			mpn_add_n(cpos, cpos, tmp, LENGTH_OF_P + 1);
		}
		else
		{
			mpn_add_n(cneg, cneg, tmp, LENGTH_OF_P + 1);
		}
	}
	if(mpn_cmp(cpos, cneg, LENGTH_OF_P + 1)>=0)
		mpn_sub_n(c_check, cpos, cneg, LENGTH_OF_P + 1);
	else
		mpn_sub_n(c_check, cneg, cpos, LENGTH_OF_P + 1);
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

/*void schoolbookn(__int128* rop, const int64_t* vect, const int64_t* matr)
SCHOOLBOOK(N)*/

_Bool barrett_equality_test(const restrict poly A, const restrict poly B)
{
	int64_t T[N];
	__int128 tmp[N];
	int64_t ctab[N]; _poly C; C.deg = N; C.t = ctab;
	// C  = A - B
	for(int i = 0; i < N; i++)
		C.t[i] = A->t[i] - B->t[i];
	
	// T <- round((C * G')/PHI)
	for(int i = 0; i < N; i++)
	{
		tmp[i] = 0;
		for(int j = 0; j < N; j++)
			tmp[i] += (__int128)C.t[j] * Gprime[j][i];
	}	
	
	for(int i = 0; i < N; i++)
		T[i] = (tmp[i] + (1ULL<<63))>>64;
	// S <- T * G
	int64_t S;
	_Bool flag = 1;
	int i = 0;
	while(i < N && flag)
	{
		S = 0;
		for(int j = 0; j < N; j++)
			S += T[j] * G[j][i];
		flag = S == C.t[i];
		i++;
	}
	return flag;
}

_Bool barrett_equality_doubletest(const __int128 * restrict A, const __int128 * restrict B)
{
	int64_t T[N];
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
	for(int i = 0; i < N; i++)
	{
		tmp[i] = 0;
		for(int j = 0; j < N; j++)
			tmp[i] += (__int128)C.t[j] * Gprime[j][i];
	}
	for(int i = 0; i < N; i++)
		T[i] = (tmp[i] + (1ULL<<63))>>64;
	
	// S <- T * G
	int64_t S;
	_Bool flag = 1;
	int i = 0;
	while(i < N && flag)
	{
		S = 0;
		for(int j = 0; j < N; j++)
			S += T[j] * G[j][i];
		flag = S == C.t[i];
		i++;
	}
	return flag;
}

void print64tab(int64_t tab[], uint16_t size, char* name)
{
	printf("%s = [", name);
	for(int i = 0; i < size-1; i++)
	{
		if(tab[i] >= 0)
			printf("%ld, ", (tab[i]));
		else
			printf("-%ld, ", -(tab[i]));
	}
	if(tab[size-1] >= 0)
		printf("%ld]\n", (tab[N-1]));
	else
		printf("-%ld]\n", -(tab[N-1]));
}

void print128tab(__int128 tab[], uint16_t size, char* name)
{
	printf("%s = [", name);
	for(int i = 0; i < size-1; i++)
	{
		if(tab[i] >= 0)
			printf("0x%lx%016lx, ", HIGH(tab[i]), LOW(tab[i]));
		else
			printf("-0x%lx%016lx, ", -HIGH(tab[i])-1, -LOW(tab[i]));
	}
	if(tab[size-1] >= 0)
		printf("0x%lx%016lx]\n", HIGH(tab[N-1]), LOW(tab[N-1]));
	else
		printf("-0x%lx%016lx]\n", -HIGH(tab[N-1])-1, -LOW(tab[N-1]));
}

void m1schoolbookn(int64_t* rop, const __int128* v, const int64_t* matr)
{
	int64_t vect[N];
	for(int i = 0; i < N; i++)
		vect[i] = v[i];
	M1SCHOOLBOOK(N)
}

#ifdef TRANSLATIONPOSSIBLE
_Bool translation_equality_halftest(const restrict poly A, const restrict poly B)
{
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = trans_T[i] + A->t[i] - B->t[i];
	
	int64_t T[N];
	for(int i = 0; i < N; i++)
	{
		T[i] = 0;
		for(int j = 0; j < N; j++)
			T[i] += C[j] * G1[j][i];
	}
	
	__int128 check = 0;
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		check = 0;
		for(int j = 0; j < N; j++)
			check += (__int128) T[j] * G[j][i];
		flag = (((__int128)C[i] + check)>>64) == 0;
		i++;
	}
	return flag;
}

_Bool translation_equality_test(const __int128 * restrict A, const __int128 * restrict B)
{
	__int128 Q[N];
	for(int i = 0; i < N; i++)
		Q[i] = A[i] - B[i];
	
	int64_t V[N];
	__int128 aux[N];
	m1toeplitz_vm(V, Q);
	mtoeplitz_vm(aux, V);
	
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = trans_T[i] + ((Q[i] + aux[i])>>64);
	
	int64_t T[N];
	
	for(int i = 0; i < N; i++)
	{
		T[i] = 0;
		for(int j = 0; j < N; j++)
			T[i] += C[j] * G1[j][i];
	}
	
	__int128 check = 0;
	_Bool flag = 1;
	uint16_t i = 0;
	while(flag && (i < N))
	{
		check = 0;
		for(int j = 0; j < N; j++)
			check += (__int128) T[j] * G[j][i];
		flag = (((__int128)C[i] + check)>>64) == 0;
		i++;
	}
	return flag;
}
#endif

_Bool extended_equality_halftest(const restrict poly A, const restrict poly B)
{
	int64_t C[N];
	for(int i = 0; i < N; i++)
		C[i] = A->t[i] - B->t[i];
	
	// hi <- (C*G1 mod E mod PHI^2)/PHI + 1
	__int128 tmp;
	int64_t S[N];
	_Bool flag = 1; // we assume equality
	_Bool earlyexit = 1;
	int i = 0;
	
	while(flag && i < N)
	{
		tmp = 0;
		for(int j = 0; j < N; j++)
			tmp += (__int128)C[j] * G2[j][i];
		S[i] = ((tmp >> 64) + 1);
		earlyexit = earlyexit && (S[i] < PHIOVER2G1) && (S[i] > -PHIOVER2G1);
		flag = flag && (tmp < STEECONST) && (tmp > -STEECONST);
		i++;
	}
	
	i = 0;
	while(!earlyexit && flag && (i < N))
	{
		tmp = 0;
		for(int j = 0; j < N; j++)
			tmp += S[j]*G[j][i];
		flag = flag && (tmp < (1ULL<<63)) && (tmp > -(1ULL<<63));
		i++;
	}
	return earlyexit && flag;
}

_Bool extended_equality_test(const __int128 * restrict A, const __int128 * restrict B)
{
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = A[i] - B[i];
	
	// hi <- (C*G1 mod E mod PHI^2)/PHI + 1
	__int128 tmp;
	int64_t S[N];
	_Bool flag = 1; // we assume equality
	_Bool earlyexit = 1;
	int i = 0;
	
	while(flag && i < N)
	{
		tmp = 0;
		for(int j = 0; j < N; j++)
			tmp += (__int128)C[j] * G2[j][i];
		S[i] = ((tmp >> 64) + 1);
		earlyexit = earlyexit && (S[i] < PHIOVER2G1) && (S[i] > -PHIOVER2G1);
		flag = flag && (tmp < DTEECONST) && (tmp > -DTEECONST);
		i++;
	}
	
	i = 0;
	while(!earlyexit && flag && (i < N))
	{
		tmp = 0;
		for(int j = 0; j < N; j++)
			tmp += S[j]*G[j][i];
		flag = flag && (tmp < (1ULL<<63)) && (tmp > -(1ULL<<63));
		i++;
	}
	return earlyexit && flag;
}

_Bool heuristic_equality_halftest(const restrict poly A, const restrict poly B)
{
	int64_t C[N];
	for(int i = 0; i < N; i++)
		C[i] = A->t[i] - B->t[i];
	__int128 Q;
	int64_t S;
	__int128 norminfQ = 0, norminfS = 0;
	_Bool flag = 1;
	int i = 0;
	while(flag && i < N)
	{
		Q = 0;
		for(int j = 0; j < N; j++)
			Q += G2[j][i] * C[j];
		S = ((Q + (1ULL<<63)) >> 64);
		norminfQ = (Q > 0? Q : -Q) > norminfQ ? (Q > 0? Q : -Q) : norminfQ;
		norminfS = (S > 0? S : -S) > norminfS ? (S > 0? S : -S) : norminfS;
		flag = (norminfQ <= STEECONST) && (norminfS < PHIOVER2G1);
		i++;
	}
	return flag;
}

_Bool heuristic_equality_test(const __int128 * restrict A, const __int128 * restrict B)
{
	__int128 C[N];
	for(int i = 0; i < N; i++)
		C[i] = A[i] - B[i];
	
	__int128 Q;
	int64_t S;
	__int128 norminfQ = 0, norminfS = 0;
	_Bool flag = 1;
	int i = 0;
	while(flag && i < N)
	{
		Q = 0;
		for(int j = 0; j < N; j++)
			Q += G2[j][i] * C[j];
		S = ((Q + (1ULL<<63)) >> 64);
		norminfQ = (Q > 0? Q : -Q) > norminfQ ? (Q > 0? Q : -Q) : norminfQ;
		norminfS = (S > 0? S : -S) > norminfS ? (S > 0? S : -S) : norminfS;
		flag = (norminfQ <= DTEECONST) && (norminfS < PHIOVER2G1);
		i++;
	}
	return flag;
}

int main(void)
{
	/*for(int i = 0; i < 2*N - 1; i++)
		M2[i] = (((__int128)M2hi[i]<<64) | (M2lo[i]));*/
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			G2[i][j] = (((__int128)G2hi[i][j]<<64) | (G2lo[i][j]));
	#ifdef TRANSLATIONPOSSIBLE
	for(int i = 0; i < N; i++)
		trans_T[i] = (((__int128)trans_Thi[i]<<64) | (trans_Tlo[i]));
	//print128tab(trans_T, N, "trans_T");
	#endif
	
	//time_t seed;
	//srand((unsigned) (time(&seed)));
	
	_poly A, B;
	int64_t atab[N], btab[N];
	A.deg = N; B.deg = N;
	A.t = atab; B.t = btab;
	
	int64_t dummy;
	_poly A2 = {.deg = N, .t = (int64_t[N]) {0} };
	_poly B2 = {.deg = N, .t = (int64_t[N]) {0} };
	
	randpoly(&A); randpoly(&B);
	for(int i = 0; i< N; i++)
	{
		A2.t[i] = A.t[i];
		B2.t[i] = B.t[i];
	}
	for(int i = 0; i < N; i++)
	{
		dummy = (rand()%7) - 3;
		for(int j = 0; j < N; j++)
			A2.t[j] += dummy*G[i][j];
	}
	for(int i = 0; i < N; i++)
	{
		dummy = (rand()%7) - 3;
		for(int j = 0; j < N; j++)
			B2.t[j] += dummy*G[i][j];
	}
	/*printf("naive\t%s", naive_equality_test(&A, &B) ? "True" : "False");
	printf("\t%s", naive_equality_test(&A, &A2) ? "True" : "False");
	printf("\t%s\n", naive_equality_test(&B, &B2) ? "True" : "False");
	printf("barrett\t%s", barrett_equality_test(&A, &B) ? "True" : "False");
	printf("\t%s", barrett_equality_test(&A, &A2) ? "True" : "False");
	printf("\t%s\n", barrett_equality_test(&B, &B2) ? "True" : "False");
	#ifdef TRANSLATIONPOSSIBLE
	printf("transla\t%s", translation_equality_halftest(&A, &B) ? "True" : "False");
	printf("\t%s", translation_equality_halftest(&A, &A2) ? "True" : "False");
	printf("\t%s\n", translation_equality_halftest(&B, &B2) ? "True" : "False");
	#endif
	printf("plantar\t%s", extended_equality_halftest(&A, &B) ? "True" : "False");
	printf("\t%s", extended_equality_halftest(&A, &A2) ? "True" : "False");
	printf("\t%s\n", extended_equality_halftest(&B, &B2) ? "True" : "False");/**/
	
#if N > 41
	int W = 6;
#else
	int W = 66;
#endif
	
	//printf("Montgomery-like cycles %ld\n", do_bench(pmns_montg_mult,W/2));
	/*printf("naive test cycles %ld\n", do_equalitybench(naive_equality_test,W/3));
	#ifdef TRANSLATIONPOSSIBLE
	printf("translat test cycles %ld\n", do_equalitybench(translation_equality_halftest,W/3));
	#else
	printf("translat test cycles -\n");
	#endif
	printf("reduced test cycles %ld\n", do_equalitybench(barrett_equality_test,W/3));
	printf("extended test cycles %ld\n", do_equalitybench(extended_equality_halftest,W/3));/**/
	
	__int128 dA[N], dB[N], dA2[N], dB2[N];
	for(int i = 0; i < N; i++)
	{
		dA[i] = randomint128() % WRHOCARRE;
		dA2[i] = dA[i];
		dB[i] = randomint128() % WRHOCARRE;
		dB2[i] = dB[i];
	}
	for(int i = 0; i < N; i++)
	{
		dummy = rand();
		for(int j = 0; j < N; j++)
			dA2[j] += (__int128)dummy*G[i][j];
	}
	for(int i = 0; i < N; i++)
	{
		dummy = rand();
		for(int j = 0; j < N; j++)
			dB2[j] += (__int128)dummy*G[i][j];
	}
	
	
	/*printf("double naive\t%s", naive_equality_doubletest(dA, dB) ? "True" : "False");
	printf("\t%s", naive_equality_doubletest(dA, dA2) ? "True" : "False");
	printf("\t%s\n", naive_equality_doubletest(dB, dB2) ? "True" : "False");
	printf("double barrett\t%s", barrett_equality_doubletest(dA, dB) ? "True" : "False");
	printf("\t%s", barrett_equality_doubletest(dA, dA2) ? "True" : "False");
	printf("\t%s\n", barrett_equality_doubletest(dB, dB2) ? "True" : "False");
	#ifdef TRANSLATIONPOSSIBLE
	printf("double transla\t%s", translation_equality_test(dA, dB) ? "True" : "False");
	printf("\t%s", translation_equality_test(dA, dA2) ? "True" : "False");
	printf("\t%s\n", translation_equality_test(dB, dB2) ? "True" : "False");
	#endif
	printf("double plantar\t%s", extended_equality_test(dA, dB) ? "True" : "False");
	printf("\t%s", extended_equality_test(dA, dA2) ? "True" : "False");
	printf("\t%s\n", extended_equality_test(dB, dB2) ? "True" : "False");/**/
	
	/*printf("naive doubletest cycles %ld\n", do_doubleequalitybench(naive_equality_doubletest,W/3));
	#ifdef TRANSLATIONPOSSIBLE
	printf("translat doubletest cycles %ld\n", do_doubleequalitybench(translation_equality_test,W/3));
	#else
	printf("translat test cycles -\n");
	#endif
	printf("reduced doubletest cycles %ld\n", do_doubleequalitybench(barrett_equality_doubletest,W/3));
	printf("extended doubletest cycles %ld\n", do_doubleequalitybench(extended_equality_test,W/3));/**/
	
	uint64_t ncycles, rcycles, ecycles, hcycles, tcycles = 0;
	char tstr[10];
	
	printf("================================================================================\n");
	
	//ncycles = 404;
	ncycles = do_equalitybench(naive_equality_test,W/3);
	#ifdef TRANSLATIONPOSSIBLE
	//tcycles = 404;
	tcycles = do_equalitybench(translation_equality_halftest,W/3);
	#endif
	//rcycles = 404;
	rcycles = do_equalitybench(barrett_equality_test,W/3);
	ecycles = do_equalitybench(extended_equality_halftest,W/3);
	hcycles = do_equalitybench(heuristic_equality_halftest,W/3);
	sprintf(tstr, tcycles ? "%ld" : "-", tcycles);
	printf("|\t%ld\t|\t%s\t|\t%ld\t|\t%ld\t|\t%ld\t|\n", ncycles, tstr, rcycles, ecycles, hcycles); 
	
	printf("================================================================================\n");
	
	ncycles = do_doubleequalitybench(naive_equality_doubletest,W/3);
	#ifdef TRANSLATIONPOSSIBLE
	tcycles = do_doubleequalitybench(translation_equality_test,W/3);
	#endif
	rcycles = do_doubleequalitybench(barrett_equality_doubletest,W/3);
	ecycles = do_doubleequalitybench(extended_equality_test,W/3);
	hcycles = do_doubleequalitybench(heuristic_equality_test,W/3);
	sprintf(tstr, tcycles ? "%ld" : "-", tcycles);
	printf("|\t%ld\t|\t%s\t|\t%ld\t|\t%ld\t|\t%ld\t|\n", ncycles, tstr, rcycles, ecycles, hcycles); 
	printf("================================================================================\n");
	return 0;
}


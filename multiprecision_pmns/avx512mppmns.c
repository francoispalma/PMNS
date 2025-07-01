#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>

#define NTEST 501
#define NSAMPLES 1001

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))

#include "avx512params.h"
#include "avx512toeplitzmacros.h"



__m512i MAVX[NBCHUNKS][2*N - 1];
__m512i M1AVX[2*N - 1];
#if N % 16 == 0
__m512i M1avxm0m1[N - 1];
__m512i M1avxm0m2[N - 1];
#endif
#if N % 32 == 0
__m512i M1avxm0m1m0m1[N/2 - 1];
__m512i M1avxm0m1m0m2[N/2 - 1];
__m512i M1avxm0m2m0m1[N/2 - 1];
__m512i M1avxm0m2m0m2[N/2 - 1];
#endif

typedef struct
{
	uint16_t deg;
	__m512i *t[NBCHUNKS];
} _mppoly, *mppoly;

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline __int128 randomint128(void)
{
	return (((__int128)randomint64())<<64) | randomint64();
}

void randpoly(mppoly P)
{
	// Function that generates a random polynomial with all coefficients < RHO.
	for(register uint16_t i = 0; i < P->deg/8; i++)
	{
		for(int16_t j = 0; j < NBCHUNKS-1; j++)
		{
			for(int k = 0; k < 8; k++)
				((int64_t*)&P->t[j][i])[k] = randomint64() % (1ULL<<(RHOver2));
		}
		for(int k = 0; k < 8; k++)
		{
			((int64_t*)&P->t[NBCHUNKS-1][i])[k] = randomint64() % (RHOhi + (RHOhi == 0));
		}
	}
}

#if N == 120
static inline void sschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	sSCHOOLBOOK(N/15)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/15)

static inline void smidtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP55TOP(N/3,sschoolbook,sschoolbook)

static inline void midtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP55TOP(N/3,schoolbook,sschoolbook)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP33TOP(N,midtoep,smidtoep)
#elif N == 80
static inline void sschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	sSCHOOLBOOK(N/10)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/10)

static inline void smidtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP55TOP(N/2,sschoolbook,sschoolbook)

static inline void midtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP55TOP(N/2,schoolbook,sschoolbook)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP22TOP(N,midtoep,smidtoep)
#elif N == 64
static inline void vschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vSCHOOLBOOK(N/8)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/8)

static inline void vmidtoep2(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/4,vschoolbook,vschoolbook)

static inline void midtoep2(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/4,vschoolbook,schoolbook)

static inline void vmidtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/2,vmidtoep2,vmidtoep2)

static inline void midtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/2,vmidtoep2,midtoep2)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N,vmidtoep,midtoep)
#elif N == 48
static inline void sschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	sSCHOOLBOOK(N/6)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/6)

static inline void smidtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP33TOP(N/2,sschoolbook,sschoolbook)

static inline void midtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP33TOP(N/2,schoolbook,sschoolbook)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP22TOP(N,midtoep,smidtoep)
#elif N == 40
static inline void sschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	sSCHOOLBOOK(N/5)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/5)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP55TOP(N,schoolbook,sschoolbook)
#elif N == 32
static inline void vschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vSCHOOLBOOK(N/4)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/4)

static inline void vmidtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/2,vschoolbook,vschoolbook)

static inline void midtoep(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N/2,vschoolbook,schoolbook)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vTOEP22TOP(N,vmidtoep,midtoep)
#elif N == 24
static inline void sschoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	sSCHOOLBOOK(N/3)

static inline void schoolbook(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N/3)

static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	TOEP33TOP(N,schoolbook,sschoolbook)
#else
// Toeplitz standard vector matrix multiplication
static inline void toeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	SCHOOLBOOK(N)
#endif

#if N == 120
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	pTOEP33TOP(N,midtoep,smidtoep)
#elif N == 80 || N == 48
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	pTOEP22TOP(N,midtoep,smidtoep)
#elif N == 64 || N == 32
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	vpTOEP22TOP(N,vmidtoep,midtoep)
#elif N == 40
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	pTOEP55TOP(N,schoolbook,sschoolbook)
#elif N == 24
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	pTOEP33TOP(N,schoolbook,sschoolbook)
#else
// Toeplitz version where the result is added to rop instead of erasing
static inline void ptoeplitz_vm(__m512i* restrict rophi, __m512i* restrict roplo, const __m512i* restrict vect, const __m512i* restrict matr)
	pSCHOOLBOOK(N)
#endif

// Toeplitz version where only the low part of the result matters
static inline void m1toep20(__m512i* restrict rop, const __m512i* restrict vect, const __m512i* restrict matr)
	M1SCHOOLBOOK(N)

static inline void m1toeplitz_vm(__m512i* restrict rop, const __m512i* restrict vect)
{
	m1toep20(rop, vect, (const __m512i*) M1AVX);
}


void pmns_montg_mult(restrict mppoly res, const restrict mppoly A, const restrict mppoly B)
{
	// Function that multiplies A by B using the Montgomery CIOS approach in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will also be in the PMNS such that if A(gamma) = a * phi mod p and
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	const __m512i PHIMASK = (__m512i) {(1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1, (1ULL<<RHOver2) - 1};
	__m512i T[N/8];
	__m512i Reslo[2*NBCHUNKS-1][N/8];
	__m512i Reshi[2*NBCHUNKS-1][N/8];
	__m512i m1aux[N/8];
	
	__m512i matr[NBCHUNKS][N/4];
	
	int64_t tmpB[NBCHUNKS][N/8][8];
	for(int i = 0; i < NBCHUNKS; i++)
		for(int j = 0; j < N/8; j++)
			_mm512_store_epi64(tmpB[i][j], B->t[i][j]);
	for(int i = 0; i < NBCHUNKS; i++)
	{
		for(int j = 0; j < N/8-1; j++)
			matr[i][j] = _mm512_set_epi64(tmpB[i][j+1][0] * LAMBDA, tmpB[i][j][7] * LAMBDA, tmpB[i][j][6] * LAMBDA, tmpB[i][j][5] * LAMBDA, tmpB[i][j][4] * LAMBDA, tmpB[i][j][3] * LAMBDA, tmpB[i][j][2] * LAMBDA, tmpB[i][j][1] * LAMBDA);
		matr[i][N/8-1] = _mm512_set_epi64(tmpB[i][0][0], tmpB[i][N/8-1][7] * LAMBDA, tmpB[i][N/8-1][6] * LAMBDA, tmpB[i][N/8-1][5] * LAMBDA, tmpB[i][N/8-1][4] * LAMBDA, tmpB[i][N/8-1][3] * LAMBDA, tmpB[i][N/8-1][2] * LAMBDA, tmpB[i][N/8-1][1] * LAMBDA);
		for(int j = 0; j < N/8-1; j++)
			matr[i][j+N/8] = _mm512_set_epi64(tmpB[i][j+1][0], tmpB[i][j][7], tmpB[i][j][6], tmpB[i][j][5], tmpB[i][j][4], tmpB[i][j][3], tmpB[i][j][2], tmpB[i][j][1]);
		matr[i][N/4 - 1] = _mm512_set_epi64(0, tmpB[i][N/8-1][7], tmpB[i][N/8-1][6], tmpB[i][N/8-1][5], tmpB[i][N/8-1][4], tmpB[i][N/8-1][3], tmpB[i][N/8-1][2], tmpB[i][N/8-1][1]);
	}
	
	// Res <- A * B[0] mod E
	for(int j = 0; j < NBCHUNKS; j++)
		toeplitz_vm((__m512i*)(Reshi[j]), (__m512i*)(Reslo[j]), (__m512i*)(A->t[j]), matr[0]);
	
	for(int i = 0; i < NBCHUNKS - 1; i++)
	{
		// T <- Res[0] * M1 mod E mod PHI
		m1toeplitz_vm((__m512i*)m1aux, Reslo[i]);
		for(int j = 0; j < N/8; j++)
		{
			T[j] = _mm512_and_si512(m1aux[j], PHIMASK);
		}
	
		// Res <- Res + T * M
		for(int j = 0; j < NBCHUNKS; j++)
			ptoeplitz_vm((__m512i*)(Reshi[i+j]), (__m512i*)(Reslo[i+j]), T, MAVX[j]);
	
		for(int j = 0; j < N/8; j++)
		{
			// Res <- Res / PHI
			Reslo[i+1][j] = _mm512_add_epi64(Reslo[i+1][j], _mm512_srai_epi64(Reslo[i][j], RHOver2));
			Reslo[i+1][j] = _mm512_add_epi64(Reslo[i+1][j], _mm512_slli_epi64(Reshi[i][j], 52-RHOver2));
			// Instead of shifting the whole result, we do it semantically
		}
	
		// Res <- Res + A * B[i+1] mod E
		for(int j = 0; j < NBCHUNKS - 1; j++)
			ptoeplitz_vm((__m512i*)(Reshi[i+1 + j]), (__m512i*)(Reslo[i+1 + j]), (__m512i*)(A->t[j]), matr[i+1]);
		toeplitz_vm((__m512i*)(Reshi[i+NBCHUNKS]), (__m512i*)(Reslo[i+NBCHUNKS]), (__m512i*)(A->t[NBCHUNKS - 1]), matr[i+1]);
	}
	
	// Last iteration is done outside the loop for small optimizations
	// T <- Res[0] * M1 mod E mod PHI
	m1toeplitz_vm((__m512i*)m1aux, Reslo[NBCHUNKS - 1]);
	
	for(int j = 0; j < N/8; j++)
	{
		T[j] = _mm512_and_si512(m1aux[j], PHIMASK);
	}
	
	// Res <- Res + T*M
	for(int j = 0; j < NBCHUNKS; j++)
		ptoeplitz_vm((__m512i*)(Reshi[NBCHUNKS - 1 + j]), (__m512i*)(Reslo[NBCHUNKS - 1 + j]), T, MAVX[j]);
	
	for(int j = 0; j < N/8; j++)
	{
		for(int k = NBCHUNKS; k < 2*NBCHUNKS-1; k++)
		{
			Reslo[k][j] = _mm512_add_epi64(Reslo[k][j], _mm512_srai_epi64(Reslo[k-1][j], RHOver2));
			Reslo[k][j] = _mm512_add_epi64(Reslo[k][j], _mm512_slli_epi64(Reshi[k-1][j], 52-RHOver2));
		}
	}
	
	// Lastly we reconstitute the result into a proper state
	for(int i = 0; i < N/8; i++)
	{
		for(int j = 0; j < NBCHUNKS - 1; j++)
			res->t[j][i] = _mm512_and_si512(Reslo[NBCHUNKS + j][i], PHIMASK);
		res->t[NBCHUNKS - 1][i] = _mm512_add_epi64(_mm512_slli_epi64(Reshi[2*NBCHUNKS-2][i], 52-RHOver2), _mm512_srai_epi64(Reslo[2*NBCHUNKS-2][i], RHOver2));
	}
}

static void p192_print(const restrict mppoly P)
{
	uint64_t tmp[NBCHUNKS];
	uint8_t counter;
	printf("[");
	for(int16_t i = 0; i < P->deg; i++)
	{
		counter = 0;
		int j = 0;
		while(j < NBCHUNKS-2)
		{
			if (counter < RHOver2)
				tmp[j] = (((int64_t*)&P->t[j][i/8])[i%8] >> counter) | (((int64_t*)&P->t[j+1][i/8])[i%8] << (RHOver2-counter)) | (((int64_t*)&P->t[j+2][i/8])[i%8] << (2*RHOver2-counter));
			else
				tmp[j] = (((int64_t*)&P->t[j+1][i/8])[i%8] >> (counter-RHOver2)) | (((int64_t*)&P->t[j+2][i/8])[i%8] << (2*RHOver2-counter));
			counter += 64 - RHOver2;
			j++;
		}
		if (counter < RHOver2)
			tmp[j] = (((int64_t*)&P->t[j][i/8])[i%8] >> counter) | (((int64_t*)&P->t[j+1][i/8])[i%8] << (RHOver2-counter));
		else
			tmp[j] = ((int64_t*)&P->t[j+1][i/8])[i%8] >> (counter-RHOver2);
		counter += 64 - RHOver2;
		j++;
		if (counter < RHOver2)
			tmp[NBCHUNKS-1] = ((int64_t*)&P->t[NBCHUNKS-1][i/8])[i%8] >> counter;
		else
			tmp[NBCHUNKS-1] = 0;
		if(((int64_t*)&P->t[NBCHUNKS-1][i/8])[i%8] >= 0)
		{
			printf("0x%lx", tmp[NBCHUNKS-1]);
			for(int16_t j = NBCHUNKS - 2; j >= 0; j--)
				printf("%016lx", tmp[j]);
		}
		else
		{
			if(NBCHUNKS > 1)
			{
			printf("-0x%lx", -tmp[NBCHUNKS-1]-1);
			for(int16_t j = NBCHUNKS - 2; j > 0; j--)
				printf("%016lx", -tmp[j] - 1);
			printf("%016lx", -tmp[0]);
			}
			else
				printf("-0x%lx", -tmp[0]);
		}
		if(i == P->deg - 1)
			printf("]\n");
		else
			printf(", ");
	}
}

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
	__asm__ __volatile__ (	"CPUID\n	"
			"RDTSC\n	"
			"mov %%edx, %0\n	"
			"mov %%eax, %1\n	"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

static inline uint64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n	"
			"mov %%edx, %0\n	"
			"mov %%eax, %1\n	"
			"CPUID\n	"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t do_bench(void (*pmns_mult)(restrict mppoly, const restrict mppoly, const restrict mppoly), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_mppoly a, b, c;
	__m512i atab[NBCHUNKS][N/8], btab[NBCHUNKS][N/8], ctab[NBCHUNKS][N/8];
	a.deg = N; b.deg = N; c.deg = N;
	for(int16_t j = 0; j < NBCHUNKS; j++)
	{
		a.t[j] = atab[j];
		b.t[j] = btab[j];
		c.t[j] = ctab[j];
	}
	mppoly A = &a, B = &b, C = &c;
	
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

int main(void)
{
	for(int j = 0; j < 2*N-1; j++)
		M1AVX[j] = _mm512_set_epi64(M1[(j+7)%(2*N-1)], M1[(j+6)%(2*N-1)], M1[(j+5)%(2*N-1)], M1[(j+4)%(2*N-1)], M1[(j+3)%(2*N-1)], M1[(j+2)%(2*N-1)], M1[(j+1)%(2*N-1)], M1[(j+0)%(2*N-1)]);
	
	#if N % 16 == 0
	for(int i = 0; i < N-1; i++)
	{
		M1avxm0m1[i] = _mm512_add_epi64(M1AVX[i+N/2],M1AVX[i+N]);
		M1avxm0m2[i] = _mm512_add_epi64(M1AVX[i+N/2],M1AVX[i]);
	}
	#endif
	#if N % 32 == 0
	for(int i = 0; i < N/2-1; i++)
	{
		M1avxm0m1m0m1[i] = _mm512_add_epi64(M1avxm0m1[i+N/4],M1avxm0m1[i+N/2]);
		M1avxm0m1m0m2[i] = _mm512_add_epi64(M1avxm0m1[i+N/4],M1avxm0m1[i]);
		M1avxm0m2m0m1[i] = _mm512_add_epi64(M1avxm0m2[i+N/4],M1avxm0m2[i+N/2]);
		M1avxm0m2m0m2[i] = _mm512_add_epi64(M1avxm0m2[i+N/4],M1avxm0m2[i]);
	}
	#endif
	
	
	for(int i = 0; i < NBCHUNKS; i++)
	{
		for(int j = 0; j < 2*N-1; j++)
			((int64_t*)&MAVX[i][j/8])[j%8] = M[i][j];
	}
	
	_mppoly A, B, C;
	__m512i atab[NBCHUNKS][N/8], btab[NBCHUNKS][N/8], ctab[NBCHUNKS][N/8];
	A.deg = N; B.deg = N; C.deg = N;
	for(int16_t j = 0; j < NBCHUNKS; j++)
	{
		A.t[j] = atab[j];
		B.t[j] = btab[j];
		C.t[j] = ctab[j];
	}
	
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	#ifndef NOBENCH
	printf("This work: %ld\n", do_bench(pmns_montg_mult, 33));
	#endif
	
	FILE *fpointer = fopen("avx512log", "w+");
	fclose(fpointer);
	fpointer = freopen("avx512log", "a+", stdout);
	
	for(int i = 0; i < 1000/(1 + 9*(NBCHUNKS*N>100)); i++)
	{
		randpoly(&A);
		randpoly(&B);
		p192_print(&A);
		p192_print(&B);
		pmns_montg_mult(&C, &A, &B);
		p192_print(&C);
	}
	
	fclose(fpointer);
	return 0;
}


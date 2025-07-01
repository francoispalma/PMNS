#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define NTEST 501
#define NSAMPLES 1001

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))

#include "toeplitzmacros.h"

#include "mpparams.h"

typedef struct
{
	uint16_t deg;
	int64_t *t[NBCHUNKS];
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
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		for(int16_t j = 0; j < NBCHUNKS-1; j++)
			P->t[j][i] = randomint64() % (1ULL<<(RHOver2));
		P->t[NBCHUNKS-1][i] = randomint64() % (RHOhi + (RHOhi == 0));
		if(randomint64() & 1)
			P->t[NBCHUNKS-1][i] *= -1;
	}
}

void pmns_montg_mult(restrict mppoly res, const restrict mppoly A, const restrict mppoly B)
{
	// Function that multiplies A by B using the Montgomery CIOS approach in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will be in the PMNS also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	int64_t T[N];
	__int128 Res[2*NBCHUNKS-1][N];
	__int128 aux[NBCHUNKS][N];
	int64_t m1aux[N];
	
	int64_t matr[NBCHUNKS][2*N-1];
	
	// We construct the Toeplitz Matrices
	for(int i = 0; i < N; i++)
		for(int j = 0; j < NBCHUNKS; j++)
			matr[j][i + N - 1] = B->t[j][i];
	for(int i = 0; i < N-1; i++)
		for(int j = 0; j < NBCHUNKS; j++)
			matr[j][i] = B->t[j][1 + i] * LAMBDA;
	
	// Res <- A * B[0] mod E
	for(int j = 0; j < NBCHUNKS; j++)
		toeplitz_vm(Res[j], A->t[j], matr[0]);
	
	for(int i = 0; i < NBCHUNKS - 1; i++)
	{
		// T <- Res[0] * M1 mod E mod PHI
		m1toeplitz_vm(m1aux, Res[i]);
		for(int j = 0; j < N; j++)
			T[j] = m1aux[j] % (1ULL<<RHOver2);
		
		// Res <- Res + T * M
		for(int j = 0; j < NBCHUNKS; j++)
			multbym[j](aux[j], (int64_t*)T);
		
		for(int j = 0; j < N; j++)
		{
			for(int k = 0; k < NBCHUNKS; k++)
				Res[i + k][j] += aux[k][j];
			// Res <- Res / PHI
			Res[i+1][j] += (Res[i][j]>>RHOver2);
			// Instead of shifting the whole result, we do it semantically
		}
		
		// Res <- Res + A * B[i+1] mod E
		for(int j = 0; j < NBCHUNKS - 1; j++)
			ptoeplitz_vm(Res[i+1 + j], A->t[j], matr[i+1]);
		toeplitz_vm(Res[i+NBCHUNKS], A->t[NBCHUNKS - 1], matr[i+1]);
	}
	
	// The last iteration is done outside the loop for small optimizations
	// T <- Res[0] * M1 mod E mod PHI
	m1toeplitz_vm(m1aux, Res[NBCHUNKS - 1]);
	for(int j = 0; j < N; j++)
		T[j] = m1aux[j] % (1ULL<<RHOver2);
	
	// Res <- Res + T*M
	for(int j = 0; j < NBCHUNKS; j++)
		multbym[j](aux[j], (int64_t*)T);
	
	for(int j = 0; j < N; j++)
	{
		Res[NBCHUNKS - 1][j] += aux[0][j];
		// We can fuse the division with the addition
		for(int k = 1; k < NBCHUNKS; k++)
			Res[NBCHUNKS - 1 + k][j] += aux[k][j] + (Res[NBCHUNKS - 2 + k][j]>>RHOver2);
	}
	
	// Lastly we reconstitute the result into a proper state
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < NBCHUNKS - 1; j++)
			res->t[j][i] = (unsigned __int128)Res[NBCHUNKS + j][i] % (1ULL<<RHOver2);
		res->t[NBCHUNKS - 1][i] = Res[2*NBCHUNKS-2][i]>>RHOver2;
	}
}

static void p192_print(const restrict mppoly P)
{
	// Printing function for multiprecision coefficient polynomials.
	uint64_t tmp[NBCHUNKS];
	uint8_t counter;
	printf("[");
	for(int16_t i = 0; i < P->deg; i++)
	{
		counter = 0;
		for(int16_t j = 0; j < NBCHUNKS-1; j++)
		{
			tmp[j] = (P->t[j][i] >> counter) | (P->t[j+1][i] << (RHOver2-counter));
			counter += 64 - RHOver2;
			counter = counter % 64;
		}
		tmp[NBCHUNKS-1] = P->t[NBCHUNKS-1][i] >> ((NBCHUNKS-1)*(64-RHOver2));
		if(P->t[NBCHUNKS-1][i] >= 0)
		{
			printf("0x%lx", tmp[NBCHUNKS-1]);
			for(int16_t j = NBCHUNKS - 2; j >= 0; j--)
				printf("%016lx", tmp[j]);
		}
		else
		{
			if(NBCHUNKS>1)
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

uint64_t do_bench(void (*pmns_mult)(restrict mppoly, const restrict mppoly, const restrict mppoly), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	_mppoly a, b, c;
	int64_t atab[NBCHUNKS][N], btab[NBCHUNKS][N], ctab[NBCHUNKS][N];
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

/****

	End of section.

*****/

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	
	#ifndef NOBENCH
	printf("This work s = %d: %ld\n", NBCHUNKS, do_bench(pmns_montg_mult,N>40 ? 3 : 33));
	#endif
	
	FILE *fpointer = fopen("log", "w+");
	fclose(fpointer);
	fpointer = freopen("log", "a+", stdout);
	
	_mppoly A, B, C;
	int64_t atab[NBCHUNKS][N], btab[NBCHUNKS][N], ctab[NBCHUNKS][N];
	A.deg = N; B.deg = N; C.deg = N;
	for(int16_t j = 0; j < NBCHUNKS; j++)
	{
		A.t[j] = atab[j];
		B.t[j] = btab[j];
		C.t[j] = ctab[j];
	}
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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "pmns256.h"
#include "pmns512.h"
#include "pmns1024.h"
#include "pmns2048.h"
#include "pmns4096.h"
#include "pmns8192.h"
#include "pmns256128.h"
#include "pmns512128.h"
#include "pmns1024128.h"
#include "pmns2048128.h"
#include "pmns4096128.h"
#include "pmns8192128.h"

#define NTEST 501
#define NSAMPLES 1001


int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

int64_t __modrho(int64_t param, const uint8_t RHO)
{
	return param & ((1ULL<<RHO) - 1);
}

void randpoly(poly P, const uint8_t RHO)
{
	// Function that generates a random polynomial with all coefficients < 2^RHO.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64(), RHO) * (1 + (rand() & 1) * -2);
}

static inline int64_t __modrhohi(int64_t param, const uint8_t RHO)
{
	// Utility function to get a high part for a random poly128.
	if (RHO <= 64)
		return 0;
	else
		return param & ((1ULL << (RHO - 64)) - 1);
}

void randpoly128(poly128 P, const uint8_t RHO)
{
	// Generates a random poly128 with appropriate, lower than rho high part.
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->lo[i] = randomint64();
		if(RHO > 64) P->hi[i] = __modrhohi(randomint64(), RHO) * (1 + (rand() & 1) * -2);
	}
}

// NTEST*NSAMPLES must be odd
// it's easier to compute median value


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

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint8_t N, const uint8_t RHO, const uint8_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(a, RHO);
		randpoly(b, RHO);
		pmns_mult(c, a, b);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(a, RHO);
		randpoly(b, RHO);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W; soak++)
				pmns_mult(c, a, b);
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
	
	free(cycles);
	free_polys(a, b, c, NULL);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_bench128(void (*pmns128_mult)(restrict poly128, const restrict poly128, const restrict poly128), const uint8_t N, const uint8_t RHO, const uint8_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	poly128 a, b, c;
	init_poly128s(N, &a, &b, &c, NULL);
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly128(a, RHO);
		randpoly128(b, RHO);
		pmns128_mult(c, a, b);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly128(a, RHO);
		randpoly128(b, RHO);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W; soak++)
				pmns128_mult(c, a, b);
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
	
	free(cycles);
	free_poly128s(a, b, c, NULL);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

int main(void)
{
	
	const uint8_t N256 = 5, N512 = 10, N1024 = 19, N2048 = 40, N4096 = 84, N8192 = 187;
	const uint8_t RHO256 = 54, RHO512 = 54, RHO1024 = 58, RHO2048 = 57, RHO4096 = 56, RHO8192 = 55;
	uint64_t cycles256, cycles512, cycles1024, cycles2048, cycles4096, cycles8192;
	const uint8_t N256128 = 3, N512128 = 5, N1024128 = 9, N2048128 = 18, N4096128 = 36, N8192128 = 72;
	const uint8_t RHO256128 = 87, RHO512128 = 105, RHO1024128 = 117, RHO2048128 = 118, RHO4096128 = 119, RHO8192128 = 120;
	uint64_t cycles256128, cycles512128, cycles1024128, cycles2048128, cycles4096128, cycles8192128;
	
	cycles256 = do_bench(pmns256_montg_mult, N256, RHO256, 10);
	cycles256128 = do_bench128(pmns256128_montg_mult, N256128, RHO256128, 10);
	printf("256-bit done\n");
	cycles512 = do_bench(pmns512_montg_mult, N512, RHO512, 10);
	cycles512128 = do_bench128(pmns512128_montg_mult, N512128, RHO512128, 10);
	printf("512-bit done\n");
	cycles1024 = do_bench(pmns1024_montg_mult, N1024, RHO1024, 10);
	cycles1024128 = do_bench128(pmns1024128_montg_mult, N1024128, RHO1024128, 10);
	printf("1024-bit done\n");
	cycles2048 = do_bench(pmns2048_montg_mult, N2048, RHO2048, 1);
	cycles2048128 = do_bench128(pmns2048128_montg_mult, N2048128, RHO2048128, 1);
	printf("2048-bit done\n");
	cycles4096 = do_bench(pmns4096_montg_mult, N4096, RHO4096, 1);
	cycles4096128 = do_bench128(pmns4096128_montg_mult, N4096128, RHO4096128, 1);
	printf("4096-bit done\n");
	cycles8192 = do_bench(pmns8192_montg_mult, N8192, RHO8192, 1);
	cycles8192128 = do_bench128(pmns8192128_montg_mult, N8192128, RHO8192128, 1);
	
	printf("\n================================================================================\n");
	printf("|  size of p   |   256   |   512   |   1024   |   2048   |   4096   |    8192  |\n");
	printf("================================================================================\n");
	printf("|    Red-64    |   %lu   |   %lu   |   %lu   |   %lu   |  %lu   |  %lu  |\n", cycles256, cycles512, cycles1024, cycles2048, cycles4096, cycles8192);
	printf("================================================================================\n");
	printf("|    Red-128   |   %lu   |   %lu   |   %lu   |  %lu   |  %lu   |  %lu  |\n", cycles256128, cycles512128, cycles1024128, cycles2048128, cycles4096128, cycles8192128);
	printf("================================================================================\n\n");
	return 0;
}


#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "pmns2048.h"
#include "pmns4096.h"
#include "pmns8192.h"

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
	__asm__ __volatile__ (	"CPUID\n\t"
			"RDTSC\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

static inline uint64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			"CPUID\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint8_t N, const uint8_t RHO)
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
			// We call the function 10 times to get an accurate measurement.
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
			pmns_mult(c, a, b);
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
	return medianTimer/NSAMPLES/10; // We divide by 10 since we measured 10 calls.
}

int main(void)
{
	const uint8_t N2048 = 40, N4096 = 84, N8192 = 187;
	const uint8_t RHO2048 = 57, RHO4096 = 56, RHO8192 = 55;
	uint64_t polycycles2048, lattcycles2048, polycycles4096, lattcycles4096, polycycles8192, lattcycles8192;
	
	printf("\nStarting.\nWARNING: The measures might be faulty if this is not run on a computer with 8 cores available.\n\n");
	polycycles2048 = do_bench(poly_pmns2048_montg_mult, N2048, RHO2048);
	lattcycles2048 = do_bench(latt_pmns2048_montg_mult, N2048, RHO2048);
	printf("2048-bit done.\n");
	polycycles4096 = do_bench(poly_pmns4096_montg_mult, N4096, RHO4096);
	lattcycles4096 = do_bench(latt_pmns4096_montg_mult, N4096, RHO4096);
	printf("4096-bit done.\n");
	polycycles8192 = do_bench(poly_pmns8192_montg_mult, N8192, RHO8192);
	lattcycles8192 = do_bench(latt_pmns8192_montg_mult, N8192, RHO8192);
	
	printf("\n================================================================================\n");
	printf("|\t\tSize\t\t|\t2048\t|\t4096\t|\t8192\t|\n");
	printf("================================================================================\n");
	printf("|\t\t\t\t\t8 cores\t\t\t\t\t|\n");
	printf("================================================================================\n");
	printf("| Toeplitz Montgomery-like\t|\t%lu\t|\t%lu\t|\t%lu\t|\n", polycycles2048, polycycles4096, polycycles8192);
	printf("================================================================================\n");
	printf("|\t\tThis work\t|\t%lu\t|\t%lu\t|\t%lu\t|\n", lattcycles2048, lattcycles4096, lattcycles8192);
	printf("================================================================================\n\n");
	return 0;
}


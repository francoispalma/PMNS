#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "pmns.h"
#include "structs.h"

#define NTEST 501
#define NSAMPLES 1001

// NTEST*NSAMPLES must be odd
// it's easier to compute median value

#ifndef RDPMC_ALLOWED

unsigned long rdpmc_instructions(void) { return 1;}

#else

inline static unsigned long rdpmc_instructions(void)
{
	unsigned a, d, c;
	
	c = (1<<30);
	__asm__ __volatile__("rdpmc" : "=a" (a), "=d" (d) : "c" (c));
	
	return ((unsigned long)a) | (((unsigned long)d) << 32);;
}

#endif

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

void do_bench(uint64_t retcycles[3], void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint8_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	poly a, b, c;
	init_polys(PDEGREE, &a, &b, &c, NULL);
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(a);
		randpoly(b);
		pmns_mult(c, a, b);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(a);
		randpoly(b);
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
	// We divide by W since we measured W calls.
	retcycles[0] = meanTimermin/NSAMPLES/W;
	retcycles[1] = meanTimermax/NSAMPLES/W;
	retcycles[2] = medianTimer/NSAMPLES/W;
}

uint64_t do_instructions(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint8_t W)
{
	uint64_t timer = 0, START, STOP, mini = (uint64_t)-1LL;
	poly a, b, c;
	init_polys(PDEGREE, &a, &b, &c, NULL);
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly(a);
		randpoly(b);
		pmns_mult(c, a, b);
	}
	
	for(int k=0;k<NSAMPLES;k++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(a);
		randpoly(b);
		for(int i=0;i<NTEST;i++)
		{
			pmns_mult(c, a, b);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			for(int soak=0; soak < W; soak++)
				pmns_mult(c, a, b);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;
		}
		timer += mini;
	}
	free_polys(a, b, c, NULL);
	return timer/NSAMPLES/W; // We divide by W since we measured W calls.
}

int main(void)
{
	//time_t seed;
	//srand((unsigned) (time(&seed)));
	
	uint64_t cycles[3];
	
	#if (PDEGREE <= 15)
		const uint8_t W = 10;
	#else
		const uint8_t W = 1;
	#endif
	
	do_bench(cycles, amns_montg_mult, W);
	
	uint64_t timer = do_instructions(amns_montg_mult, W);
	
	printf("(%ld, %ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2], timer);
	return 0;
}


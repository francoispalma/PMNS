#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <time.h>
#include <inttypes.h>

#include "structs.h"
#include "pmns256.h"
#include "pmns512.h"
#include "pmns1024.h"

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

unsigned long long *quartiles(uint64_t *tab, uint64_t size)
{
	unsigned long long *result ;
	uint64_t aux ;
	
	result = malloc(3*sizeof(unsigned long long));
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


int main(void)
{
	uint64_t *cycles1 ;
	unsigned long long timermin1 , timermax1, meanTimer1min =0,	medianTimer1 = 0,
	meanTimer1max = 0, t1,t2, diff_t;
	
	unsigned long long *statTimer1 ;
	
	uint64_t polycycles256, lattcycles256, polycycles512, lattcycles512, polycycles1024, lattcycles1024;
	const uint8_t N256 = 5, N512 = 10, N1024 = 19;
	const uint8_t RHO256 = 54, RHO512 = 54, RHO1024 = 58;
	uint8_t N, RHO;
	
	poly a, b, c;
	
	void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly);
	
	N = N256; RHO = RHO256; 	
	
	init_polys(N, &a, &b, &c, NULL);
	
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = poly_pmns256_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	polycycles256 = medianTimer1/NSAMPLES/10;
	
	free(cycles1);
	
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = latt_pmns256_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	lattcycles256 = medianTimer1/NSAMPLES/10;
	
	free_polys(a, b, c, NULL);
	free(cycles1);
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	
	
	/*************************************************
	512 bits
	*************************************************/
	
	N = N512; RHO = RHO512;
	
	init_polys(N, &a, &b, &c, NULL);
	
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = poly_pmns512_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	polycycles512 = medianTimer1/NSAMPLES/10;
	
	free(cycles1);
	
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = latt_pmns512_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	lattcycles512 = medianTimer1/NSAMPLES/10;
	
	free_polys(a, b, c, NULL);
	free(cycles1);
	
	/*************************************************
	1024 bits
	*************************************************/
	
	N = N1024; RHO = RHO1024;
	
	init_polys(N, &a, &b, &c, NULL);
	
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = poly_pmns1024_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	polycycles1024 = medianTimer1/NSAMPLES/10;
	
	free(cycles1);
	
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	pmns_mult = latt_pmns1024_montg_mult;
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
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
		memset(cycles1,0,NTEST*sizeof(uint64_t));
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
			if(timermin1>diff_t) timermin1 = diff_t;
			else if(timermax1 < diff_t) timermax1 = diff_t;
			cycles1[j]=diff_t;
		}
		meanTimer1min += timermin1;
		meanTimer1max += timermax1;
		statTimer1	 = quartiles(cycles1,NTEST);
		medianTimer1 += statTimer1[1];
		free(statTimer1);
	}
	
	// We divide by 10 since we measured 10 calls.
	lattcycles1024 = medianTimer1/NSAMPLES/10;
	
	free_polys(a, b, c, NULL);
	free(cycles1);
	meanTimer1min = 0; medianTimer1 = 0; meanTimer1max = 0;
	printf("\n=========================================================================\n");
	printf("|\tAlg.\\size of p\t|\t256\t|\t512\t|\t1024\t|\n");
	printf("=========================================================================\n");
	printf("|\tMontgomery\t|\t%lu\t|\t%lu\t|\t%lu\t|\n", polycycles256, polycycles512, polycycles1024);
	printf("|\tThis work\t|\t%lu\t|\t%lu\t|\t%lu\t|\n", lattcycles256, lattcycles512, lattcycles1024);
	printf("=========================================================================\n\n");
	return 0;
}


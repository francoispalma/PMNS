#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <immintrin.h>
#include <time.h>
#include <inttypes.h>

#include "structs.h"
#include "params.h"

#define NTEST 51
#define NSAMPLES 101

/*inline static unsigned long rdpmc_instructions(void)*/
/*{*/
/*   unsigned a, d, c;*/

/*   c = (1<<30);*/
/*   __asm__ __volatile__("rdpmc" : "=a" (a), "=d" (d) : "c" (c));*/

/*   return ((unsigned long)a) | (((unsigned long)d) << 32);;*/
/*}*/

unsigned long rdpmc_instructions(void) { return 1;}

/*extern void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
extern void amns128_montg_mult_pre(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);
extern void amns128_montg_mult_hyb(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);*/
extern void amns_ltr_sqandmult(restrict poly res, const restrict poly base,
	const restrict mpnum exponent);
void amns_montg_ladder(restrict poly res, const restrict poly base,
	const restrict mpnum exponent);
void randpoly(poly);

// NTEST*NSAMPLES must be odd
// it's easier to compute median value


/**** Measurements procedures according to INTEL white paper

  "How to benchmark code execution times on INTEL IA-32 and IA-64"

 *****/

void quicksort(uint64_t *tab, int n);
unsigned long long *quartiles(uint64_t *tab, uint64_t size);

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
	result[1]  = tab[(size+1)/2 - 1];
	// Q3
	aux = (3*size) >> 2;
	if ((3*size) % 4) aux++;
	result[2]  = tab[aux - 1];

	return result;
}

inline static uint64_t cpucyclesStart (void) {

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

inline static uint64_t cpucyclesStop (void) {

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


static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

void randmpnumber(mpnum P)
{
	for(register uint16_t i = 0; i < P->deg - 1; i++)
		P->t[i] = randomint64();
	P->t[P->deg - 1] = __modrho(randomint64());
	P->sign = 1;
}

int main(int argc, char** argv)
{
  uint64_t *cycles1 ;
  unsigned long long timermin1 , timermax1, meanTimer1min =0,  medianTimer1 = 0,
  meanTimer1max = 0, t1,t2, diff_t;

  unsigned long long *statTimer1 ;

	unsigned long long int timer = 0;
	uint64_t mini = (uint64_t)-1L;
	unsigned long long int START, STOP;

	poly a, c, soak1, soak2;
	mpnum b;
	init_polys(N, &a, &b, &c, &soak1, &soak2, NULL);
	init_mpnum(N, &b);
	
	void (*amns_exp)(restrict poly, const restrict poly,
		const restrict mpnum);
	if (argc > 1 && (strncmp(argv[1], "sqm", 3) == 0))
		amns_exp = amns_ltr_sqandmult;
	else
		amns_exp = amns_montg_ladder;

  cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));

  for(int i=0;i<NTEST;i++)
  {
    // ici tirage de donnees aleatoires
    // et execution de la fonction a mesurer
    // pour chauffer les caches
		randpoly(a);
		randmpnumber(b);
		amns_exp(c, a, b);
  }

  for(int i=0;i<NSAMPLES;i++)
	{
		// generer ici un jeu de parametres aleatoire pour la
		// fonction a mesurer
		randpoly(a);
		randmpnumber(b);
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
        memset(cycles1,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
            // appel de la fonction a mesurer
			amns_exp(c, a, b);
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
        statTimer1   = quartiles(cycles1,NTEST);
        medianTimer1 += statTimer1[1];
        free(statTimer1);
	}
	
	timer=0;
	
	for(int k=0; k<NSAMPLES;k++)
	{
		for(int i=0;i<NTEST;i++)
		{
			amns_exp(c, a, b);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			amns_exp(c, a, b);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;
		}
		timer += mini;
	}
	
	printf("(%lld, %lld, %lld, %lld)\n", meanTimer1min/NSAMPLES, meanTimer1max/NSAMPLES, medianTimer1/NSAMPLES, timer/NSAMPLES);
	free_polys(a, c, soak1, soak2, NULL);
	free_mpnum(b);
	free(cycles1);
	return 0;

}

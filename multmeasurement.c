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

#define NTEST 511
#define NSAMPLES 1001

/*inline static unsigned long rdpmc_instructions(void)*/
/*{*/
/*   unsigned a, d, c;*/

/*   c = (1<<30);*/
/*   __asm__ __volatile__("rdpmc" : "=a" (a), "=d" (d) : "c" (c));*/

/*   return ((unsigned long)a) | (((unsigned long)d) << 32);;*/
/*}*/

unsigned long rdpmc_instructions(void) { return 1;}

extern void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B);
extern void m1_mns_mod_mult_ext_red(int64_t* restrict R,
	const restrict poly A);
void randpoly(poly);

void mns_multab_andmultm1(int64_t* restrict V, __int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	poly tmp;
	init_poly(N, &tmp);
	
	mns_mod_mult_ext_red(R, A, B);
	
	for(register int16_t i = 0; i < N; i++)
	{
		V[i] = R[i];
		tmp->t[i] = R[i];
		R[i] = 0;
	}
	
	m1_mns_mod_mult_ext_red(V, tmp);
	
	free_poly(tmp);
}

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


int main(int argc, char** argv)
{
  uint64_t *cycles1 ;
  unsigned long long timermin1 , timermax1, meanTimer1min =0,  medianTimer1 = 0,
  meanTimer1max = 0, t1,t2, diff_t;

  unsigned long long *statTimer1 ;

	unsigned long long int timer = 0;
	uint64_t mini = (uint64_t)-1L;
	unsigned long long int START, STOP;
	
	poly a, b, soak1, soak2;
	init_polys(N, &a, &b, &soak1, &soak2, NULL);
	__int128 c[N];
	int64_t V[N];
	
	void (*mult_func)(__int128* restrict, const restrict poly, const restrict poly);
	if (argc > 1 && (strncmp(argv[1], "pre", 3) == 0))
		mult_func = mns_mod_mult_ext_red_pre;
	else
		mult_func = mns_mod_mult_ext_red;

  cycles1 = (uint64_t *)calloc(NTEST,sizeof(uint64_t));

  for(int i=0;i<NTEST;i++)
  {
    // ici tirage de donnees aleatoires
    // et execution de la fonction a mesurer
    // pour chauffer les caches
		randpoly(a);
		randpoly(b);
		mult_func(c, a, b);
		//mns_multab_andmultm1(V, c, a, b);
  }

  for(int i=0;i<NSAMPLES;i++)
	{
		// generer ici un jeu de parametres aleatoire pour la
		// fonction a mesurer
		randpoly(a);
		randpoly(b);
		for(int j=0;j<NTEST;j++)
		{
			mult_func(c, a, b);
			//mns_multab_andmultm1(V, c, a, b);
		}
		timermin1 = (unsigned long long int)0x1<<63;
		timermax1 = 0;
        memset(cycles1,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
            // appel de la fonction a mesurer
			mult_func(c, a, b);
			//mns_multab_andmultm1(V, c, a, b);
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
			mult_func(c, a, b);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mult_func(c, a, b);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;
		}
		timer += mini;
	}
	
	printf("(%lld, %lld, %lld, %lld)\n", meanTimer1min/NSAMPLES, meanTimer1max/NSAMPLES, medianTimer1/NSAMPLES, timer/NSAMPLES);
	free_polys(a, b, soak1, soak2, NULL);
	free(cycles1);
	return 0;
}


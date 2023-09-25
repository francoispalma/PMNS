#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <omp.h>

#define NTEST 501
#define NSAMPLES 1001


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

void mpnmulnwrapper(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b)
{
	mpn_mul_n(c, a, b, 128);
}

void lololohihilohihi(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b)
{
	mp_limb_t clolo[4][128];
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i < 4; i++)
		mpn_mul_n(clolo[i], a + 64 *((i & 2) > 0), b + 64 * (i & 1), 64);
	c[128] += mpn_add_n(c, c, clolo[0], 128);
	c[192] += mpn_add_n(c + 64, c + 64, clolo[1], 128);
	c[192] += mpn_add_n(c + 64, c + 64, clolo[2], 128);
	mpn_add_n(c + 128, c + 128, clolo[3], 128);
}

static inline uint64_t gmpbench(mpz_t A, mpz_t B, mp_limb_t *c_limbs, gmp_randstate_t r, uint8_t W, void (*gmpmul)(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b))
{
	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1, t2, diff_t, *statTimer;
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	mp_limb_t *a_limbs, *b_limbs;
	int nb_limbs = 128;
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		mpz_urandomb(A, r, 8192);
		mpz_urandomb(B, r, 8192);
		
		a_limbs = mpz_limbs_modify (A, nb_limbs);
		b_limbs = mpz_limbs_modify (B, nb_limbs);
		gmpmul(c_limbs, a_limbs, b_limbs);
	}
	
	{
		mp_limb_t tempc[256], sub[256];
		for(int i = 0; i < 256; i++)
		{
			c_limbs[i] = 0;
			tempc[i] = 0;
		}
		mpnmulnwrapper(c_limbs, a_limbs, b_limbs);
		lololohihilohihi(tempc, a_limbs, b_limbs);
		//gmp_printf("%Nx\n\n%Nx\n\n", c_limbs, 256, tempc, 256);
		mpn_sub_n(sub, tempc, c_limbs, 256);
		gmp_printf("0x%Nx\n", sub, 256);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		mpz_urandomb(A, r, 8192);
		mpz_urandomb(B, r, 8192);
		
		a_limbs = mpz_limbs_modify (A, nb_limbs);
		b_limbs = mpz_limbs_modify (B, nb_limbs);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W; soak++)
				gmpmul(c_limbs, a_limbs, b_limbs);
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
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}



void do_benchgmp(uint64_t retcycles[3])
{
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	int nb_limbs = 128;
	
	mpz_t A, B;
	mpz_inits(A, B, NULL);
	
	mp_limb_t *c_limbs;
	
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	
	retcycles[1] = gmpbench(A, B, c_limbs, r, 5, lololohihilohihi);
	retcycles[0] = gmpbench(A, B, c_limbs, r, 5, mpnmulnwrapper);
	
	mpz_clears(A, B, NULL);
	free(c_limbs);
	gmp_randclear(r);
}


int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cyclesGMP8192[3];
	
	do_benchgmp(cyclesGMP8192);
	
	printf("=======================================\n");
	printf("|  Low level  |   %lu   |   %lu   |\n", cyclesGMP8192[0], cyclesGMP8192[1]);
	printf("=======================================\n");
	return 0;
}

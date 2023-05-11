#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "pmns1024.h"
#include "pmns1664.h"
#include "pmns2048.h"
#include <gmp.h>
#include "gmp_stuff.c"
#include "utilitymp.h"

#define NTEST 501
#define NSAMPLES 1001

#define UNUSED(X) (void)(X)

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

static inline void gmp_lowlevel_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *c_limbs, mp_limb_t *q_limbs, int nb_limbs)
{
	mpn_mul_n(c_limbs, a_limbs, b_limbs, nb_limbs); // compute: z = y*x
	mpn_tdiv_qr(q_limbs, a_limbs, 0, c_limbs, (nb_limbs*2), p_limbs, nb_limbs); // compute: y = z%p
}

static inline void gmp_montgomery_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *mip_limbs, mp_limb_t *soak, int nb_limbs)
{
	UNUSED(soak);
	mpn_mont_mul_red_n(a_limbs, a_limbs, b_limbs, p_limbs, mip_limbs, nb_limbs);
}

static inline void gmp_montgomeryCIOS_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *mip0, mp_limb_t *soak, int nb_limbs)
{
	UNUSED(soak);
	mpn_mont_mul_red_1(a_limbs, a_limbs, b_limbs, p_limbs, *mip0, nb_limbs);
}

static inline uint64_t gmpbench(mpz_t A, mpz_t B, mpz_t modul_p, gmp_randstate_t r, mp_limb_t *c_limbs, mp_limb_t *q_limbs, uint8_t W, void (*gmp_wrapper)(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *c_limbs, mp_limb_t *q_limbs, int nb_limbs))
{
	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1, t2, diff_t, *statTimer;
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	mp_limb_t *p_limbs, *a_limbs, *b_limbs;
	int nb_limbs = mpz_size(modul_p);
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		mpz_urandomm(A, r, modul_p);
		mpz_urandomm(B, r, modul_p);
		
		p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
		a_limbs = mpz_limbs_modify (A, nb_limbs);
		b_limbs = mpz_limbs_modify (B, nb_limbs);
		gmp_wrapper(a_limbs, b_limbs, p_limbs, c_limbs, q_limbs, nb_limbs);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		mpz_urandomm(A, r, modul_p);
		mpz_urandomm(B, r, modul_p);
	
		p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
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
				gmp_wrapper(a_limbs, b_limbs, p_limbs, c_limbs, q_limbs, nb_limbs);
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



void do_benchgmp(uint64_t retcycles[3], mpnum Pp, const uint8_t W)
{
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	mpz_t modul_p;
	mpz_init(modul_p);
	
	modul_p->_mp_alloc = Pp->deg;
	modul_p->_mp_size = Pp->deg * Pp->sign;
	modul_p->_mp_d = Pp->t;
	
	int nb_limbs = mpz_size(modul_p);
	
	mpz_t A, B;
	mpz_inits(A, B, NULL);
	
	mp_limb_t *p_limbs, *c_limbs, *q_limbs, *mip_limbs;
	
	q_limbs = (mp_limb_t*) calloc ((nb_limbs+1), sizeof(mp_limb_t));
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	
	retcycles[0] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_lowlevel_wrapper);
	
	p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	mpn_binvert(mip_limbs, p_limbs, nb_limbs, c_limbs);
	
	retcycles[1] = gmpbench(A, B, modul_p, r, mip_limbs, q_limbs, W, gmp_montgomery_wrapper);
	
	mp_limb_t mip0;
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	retcycles[2] = gmpbench(A, B, modul_p, r, &mip0, q_limbs, W, gmp_montgomeryCIOS_wrapper);
	
	modul_p->_mp_alloc = 0;
	modul_p->_mp_size = 0;
	modul_p->_mp_d = NULL;
	
	mpz_clears(modul_p, A, B, NULL);
	free(c_limbs);
	free(q_limbs);
	free(mip_limbs);
	gmp_randclear(r);
}

void horner_modulo(mpnum* res, poly P, mpnum gamma, mpnum prime)
{
	// Function that evaluates the polynomial P in gamma modulo the prime.
	// Uses the horner algorithm.
	
	mpnum tmp, soak;
	init_mpnums(1, &tmp, &soak, NULL);
	uint16_t N = P->deg;
	
	for(int i = 0; i < N - 1; i++)
	{
		tmp->sign = 1 - 2*(P->t[N - 1 - i] < 0);
		tmp->t[0] = P->t[N - 1 - i] * tmp->sign;
		mp_add(&soak, tmp, *res);
		mp_mult(res, soak, gamma);
		mp_mod(&soak, *res, prime);
		mp_copy(res, soak);
	}
	tmp->sign = 1 - 2*(P->t[0] < 0);
	tmp->t[0] = P->t[0] * tmp->sign;
	mp_add(&soak, *res, tmp);
	mp_mod(res, soak, prime);
	free_mpnums(tmp, soak, NULL);
}

void check_pmns_vs_gmp(mpnum __P__, uint16_t N, mpnum gamma, uint8_t RHO, void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly))
{
	// Function to check if our computed values are correct or not.
	
	// PHI = 2^64
	_mpnum PHI = { .deg = 2, .sign = 1, .t = (uint64_t[]) {0, 1} };
	
	poly A, B, C;
	mpz_t gA, gB, gC, gC_check, gmP;
	mpnum mpA, mpB, mpC, tmp;
	mpz_inits(gA, gB, gC, gC_check, gmP, NULL);
	init_polys(N, &A, &B, &C, NULL);
	init_mpnums(__P__->deg, &mpA, &mpB, &mpC, &tmp, NULL);
	
	// We set a mpz_t with the value of P
	gmP->_mp_size = __P__->deg * __P__->sign;
	gmP->_mp_alloc = __P__->deg;
	gmP->_mp_d = __P__->t;
	
	// We generate two polynomials in our PMNS
	randpoly(A, RHO);
	randpoly(B, RHO);
	
	// We evaluate them in gamma to see which integer they represent
	horner_modulo(&mpA, A, gamma, __P__);
	horner_modulo(&mpB, B, gamma, __P__);
	
	// We transfer the values to mpz_t variables for later
	gA->_mp_size = mpA->deg * mpA->sign;
	gA->_mp_alloc = mpA->deg;
	gA->_mp_d = mpA->t;
	
	gB->_mp_size = mpB->deg * mpB->sign;
	gB->_mp_alloc = mpB->deg;
	gB->_mp_d = mpB->t;
	
	// We compute C = AB times PHI^-1 mod P
	pmns_mult(C, A, B);
	horner_modulo(&mpC, C, gamma, __P__);
	
	// We multiply by PHI mod P to get C = AB mod P
	mp_mult(&tmp, mpC, &PHI); 
	mp_mod(&mpC, tmp, __P__);
	
	// we compute C = AB mod P in GMP to check
	mpz_mul(gC, gA, gB); 
	mpz_mod(gC, gC, gmP); 
	
	if(0) // change this to 1 if you want to print the values
	{
		gmp_printf("################ A ###############\n0x%Zx\n", gA);
		mp_print(mpA);
		gmp_printf("################ B ###############\n0x%Zx\n", gB);
		mp_print(mpB);
		gmp_printf("################ C ###############\n0x%Zx\n", gC);
		mp_print(mpC);
	}
	
	// We transfer our computed C = AB mod P to check using mpz_cmp later
	gC_check->_mp_size = mpC->deg * mpC->sign;
	gC_check->_mp_alloc = mpC->deg;
	gC_check->_mp_d = mpC->t;
	
	printf("A times B GMP equals A times B PMNS: %s\n\n", 
		mpz_cmp(gC, gC_check) == 0 ? "True" : "False");
	
	free_mpnums(mpA, mpB, mpC, tmp, NULL);
	free_polys(A, B, C, NULL);
	
	// We make sure to put them back at what they were to avoid double frees
	gA->_mp_size = 0;
	gA->_mp_alloc = 0;
	gA->_mp_d = NULL;
	gB->_mp_size = 0;
	gB->_mp_alloc = 0;
	gB->_mp_d = NULL;
	gC_check->_mp_size = 0;
	gC_check->_mp_alloc = 0;
	gC_check->_mp_d = NULL;
	gmP->_mp_size = 0;
	gmP->_mp_alloc = 0;
	gmP->_mp_d = NULL;
	
	mpz_clears(gA, gB, gC, gC_check, gmP, NULL);
}

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles1024, cycles1664, cycles2048;
	uint64_t cyclesGMP1024[3], cyclesGMP1664[3], cyclesGMP2048[3];
	
	printf("\nFirst checking A times B mod P for each size of PMNS\n");
	
	printf("Checking 1024 bits:\n");
	check_pmns_vs_gmp(__P1024__, N1024, Gamma1024, RHO1024, pmns1024_montg_mult);
	
	printf("Checking 1664 bits:\n");
	check_pmns_vs_gmp(__P1664__, N1664, Gamma1664, RHO1664, pmns1664_montg_mult);
	
	printf("Checking 2048 bits:\n");
	check_pmns_vs_gmp(__P2048__, N2048, Gamma2048, RHO2048, pmns2048_montg_mult);
	
	printf("\nNext measuring execution times for each prime size\n\n");
	
	cycles1024 = do_bench(pmns1024_montg_mult, N1024, RHO1024, 10);
	do_benchgmp(cyclesGMP1024, __P1024__, 10);
	printf("1024-bit done\n");
	cycles1664 = do_bench(pmns1664_montg_mult, N1664, RHO1664, 1);
	do_benchgmp(cyclesGMP1664, __P1664__, 1);
	printf("1664-bit done\n");
	cycles2048 = do_bench(pmns2048_montg_mult, N2048, RHO2048, 1);
	do_benchgmp(cyclesGMP2048, __P2048__, 1);
	
	printf("\n=========================================================\n");
	printf("|         Size          |   1024   |   1664   |   2048  |\n");
	printf("=========================================================\n");
	printf("|       Low level       |   %lu   |   %lu   |  %lu   |\n", cyclesGMP1024[0], cyclesGMP1664[0], cyclesGMP2048[0]);
	printf("=========================================================\n");
	printf("|    Classical Mont.    |   %lu   |   %lu   |  %lu   |\n", cyclesGMP1024[1], cyclesGMP1664[1], cyclesGMP2048[1]);
	printf("=========================================================\n");
	printf("|       mont. CIOS      |   %lu   |   %lu   |  %lu   |\n", cyclesGMP1024[2], cyclesGMP1664[2], cyclesGMP2048[2]);
	printf("=========================================================\n");
	printf("| Toeplitz (this work)  |   %lu   |   %lu   |  %lu   |\n", cycles1024, cycles1664, cycles2048);
	printf("=========================================================\n");
	return 0;
}


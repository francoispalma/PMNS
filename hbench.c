#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "hparams.h"
#include "hpmns.h"

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

#define BERNLO(X) (((int64_t)X) & 4503599627370495)
#define BERNHI(X) ((int64_t)(X>>51))

#define NTEST 501
#define NSAMPLES 1001

#define KPRIMEC 189

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

void randpoly(poly P)
{
	// Function that generates a random polynomial with all coefficients < 2^RHO.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64()) * (1 + (rand() & 1) * -2);
}

void randlimbs(mp_limb_t g[])
{
	for(int i = 0; i < GMPLIMB; i++)
		g[i] = randomint64();
}

void randbernlimbs(mp_limb_t g[])
{
	for(int i = 0; i < GMPLIMB + 1; i++)
		g[i] = randomint64() & 4503599627370495;
}

void randmersennelimbs(mp_limb_t g[])
{
	for(int i = 0; i < 9; i++)
		g[i] = randomint64();
	g[8] = g[8] % 512;
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

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	
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
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				pmns_mult(c, a, b);
				pmns_mult(b, c, a);
				pmns_mult(a, b, c);
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
	
	free(cycles);
	free_polys(a, b, c, NULL);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_gmpbench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[]), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	mp_limb_t a[GMPLIMB], b[GMPLIMB], c[GMPLIMB];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randlimbs(a);
		randlimbs(b);
		gmp_mult(c, a, b);
	}
	
	gmp_printf("0x%Nx\n0x%Nx\n0x%Nx\n", a, GMPLIMB, b, GMPLIMB, c, GMPLIMB);
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randlimbs(a);
		randlimbs(b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				gmp_mult(c, a, b);
				gmp_mult(b, c, a);
				gmp_mult(a, b, c);
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
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

mp_limb_t mtmp[18];
mp_limb_t shift[9];
uint64_t do_mersennebench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[]), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	mp_limb_t a[9], b[9], c[9];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randmersennelimbs(a);
		randmersennelimbs(b);
		gmp_mult(c, a, b);
	}
	
	gmp_printf("0x%Nx\n0x%Nx\n0x%Nx\n", a, 9, b, 9, c, 9);
	gmp_printf("mtmp = 0x%Nx\n", mtmp, 18);
	gmp_printf("shift = 0x%Nx\n", shift, 9);
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randmersennelimbs(a);
		randmersennelimbs(b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W; soak++)
				gmp_mult(c, a, b);
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

void bernprint(mp_limb_t* P)
{
	printf("0x%lx", BERNLO(P[GMPLIMB]));
	for(int16_t i = GMPLIMB - 1; i > -1; i--)
		printf("%013lx", BERNLO(P[i]));
	printf("%013lx\n", BERNLO(P[0]));
}

uint64_t do_bernbench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[]), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	mp_limb_t a[GMPLIMB + 1], b[GMPLIMB + 1], c[GMPLIMB + 1];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randbernlimbs(a);
		randbernlimbs(b);
		gmp_mult(c, a, b);
	}
	
	bernprint(a);
	bernprint(b);
	bernprint(c);
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randbernlimbs(a);
		randbernlimbs(b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W; soak++)
				gmp_mult(c, a, b);
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

mp_limb_t tmp[GMPLIMB * 2];
mp_limb_t hihi[GMPLIMB], lolo[GMPLIMB], amid[GMPLIMB/2], bmid[GMPLIMB/2], mid[GMPLIMB];
void gmpmulmod2k(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	mpn_mul_n(tmp, a, b, GMPLIMB);
	/*mpn_mul_n(tmp, a, b, GMPLIMB/2);
	mpn_mul_n(tmp + GMPLIMB, a + GMPLIMB/2, b + GMPLIMB/2, GMPLIMB/2);
	mpn_add_n(amid, a, a + GMPLIMB/2, GMPLIMB/2);
	mpn_add_n(bmid, b, b + GMPLIMB/2, GMPLIMB/2);
	mpn_mul_n(mid, amid, bmid, GMPLIMB/2);
	mpn_sub_n(mid, mid, tmp + GMPLIMB, GMPLIMB);
	mpn_sub_n(mid, mid, tmp, GMPLIMB);
	mpn_add_n(tmp + GMPLIMB/2, tmp + GMPLIMB/2, mid, GMPLIMB);*/
	mpn_add_1(tmp, tmp, 1, mpn_mul_1(c, tmp + GMPLIMB, GMPLIMB, KPRIMEC) * KPRIMEC);
	mpn_add_n(c, c, tmp, GMPLIMB);
	//mpn_add_n(c, tmp, tmp + GMPLIMB, GMPLIMB);
}

__int128 tmpp[GMPLIMB + 1];
void alamano(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	__int128 aux, extra;
	tmpp[0] = 0;
	for(int i = 0; i < GMPLIMB - 1; i++)
	{
		aux = (__int128) a[i + 1] * b[GMPLIMB - 1];
		tmpp[i] += LOW(aux);
		tmpp[i + 1] = HIGH(aux);
		for(int j = 2; j < GMPLIMB - i; j++)
		{
			aux = (__int128) a[i + j] * b[GMPLIMB - j];
			tmpp[i] += LOW(aux);
			tmpp[i + 1] += HIGH(aux);
		}
		tmpp[i] *= KPRIMEC;
	}
	tmpp[GMPLIMB - 1] *= KPRIMEC;
	
	for(int i = 0; i < GMPLIMB - 1; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			aux = (__int128) a[j] * b[i - j];
			tmpp[i] += LOW(aux);
			tmpp[i + 1] += HIGH(aux);
		}
	}
	aux = (__int128) a[0] * b[GMPLIMB - 1];
	extra = HIGH(aux);
	tmpp[GMPLIMB - 1] += LOW(aux);
	for(int j = 1; j < GMPLIMB; j++)
	{
		aux = (__int128) a[j] * b[GMPLIMB - 1 - j];
		tmpp[GMPLIMB - 1] += LOW(aux);
		extra += HIGH(aux);
	}
	tmpp[0] += KPRIMEC * (extra + HIGH(tmpp[GMPLIMB - 1]));
	
	c[0] = tmpp[0];
	for(int i = 1; i < GMPLIMB; i++)
	{
		tmpp[i] += HIGH(tmpp[i - 1]);
		c[i] = tmpp[i];
	}
}



void alabern(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	for(int i = 0; i < GMPLIMB; i++)
	{
		tmpp[i] = (__int128) a[i + 1] * b[GMPLIMB];
		for(int j = 2; j < GMPLIMB + 1 - i; j++)
		{
			tmpp[i] += (__int128) a[i + j] * b[GMPLIMB + 1 - j];
		}
		tmpp[i] *= KPRIMEC;
	}
	tmpp[GMPLIMB] = 0;
	
	tmpp[0] += (__int128) a[0] * b[0];
	tmpp[1] += BERNHI(tmpp[0]);
	for(int i = 1; i < GMPLIMB; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			tmpp[i] += (__int128) a[j] * b[i - j];
		}
		tmpp[i + 1] += BERNHI(tmpp[i]);
		c[i] = BERNLO(tmpp[i]);
	}
	for(int j = 0; j < GMPLIMB + 1; j++)
	{
		tmpp[GMPLIMB] += (__int128) a[j] * b[GMPLIMB - j];
	}
	c[GMPLIMB] = BERNLO(tmpp[GMPLIMB]);
	c[0] = BERNLO(tmpp[0]) + BERNHI(tmpp[GMPLIMB]) * KPRIMEC;
	
}

void mersenne521(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	mpn_mul_n(mtmp, a, b, 9);
	mpn_rshift(shift, mtmp + 8, 9, 9);
	mpn_add_n(c, shift, mtmp, 9);
	c[8] &= 511;
}

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles;
	
	cycles = do_bench(pmns_montg_mult, 201);
	printf("pmns %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 201);
	printf("gmpmulmod2k %ld\n", cycles);
	
	/*cycles = do_gmpbench(alamano, 200);
	printf("alamano %ld\n", cycles);
	
	cycles = do_bernbench(alabern, 200);
	printf("bern %ld\n", cycles);*/
	
	/*cycles = do_mersennebench(mersenne521, 200);
	printf("%ld\n", cycles);*/
	
	return 0;
}


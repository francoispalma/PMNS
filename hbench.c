#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "hparams.h"
#include "hpmns.h"
#include "eccoptimizedcode.h"

#define LOW(X) ((uint64_t)X)
#define LO(X) ((int64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

#define BERNLO(X) (((int64_t)X) & 4503599627370495)
#define BERNHI(X) ((int64_t)(X>>51))

#define NTEST 501
#define NSAMPLES 1001

/*#define MPLIMB 9*/
/*#define KPRIMEC 17*/
/*#define LASTLIMBMASK 0x8000000000000000*/

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

/*void randlimbs(mp_limb_t g[])*/
/*{*/
/*	for(int i = 0; i < MPLIMB; i++)*/
/*		g[i] = randomint64();*/
/*}*/

void randhlimbs(mp_limb_t g[], const uint64_t limbs)
{
	for(uint64_t i = 0; i < limbs; i++)
		g[i] = randomint64();
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

uint64_t do_gmpbench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMBMASK, const mp_limb_t MPLIMB), const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMBMASK, const mp_limb_t MPLIMB, const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	mp_limb_t a[MPLIMB], b[MPLIMB], c[MPLIMB];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randhlimbs(a, MPLIMB);
		randhlimbs(b, MPLIMB);
		gmp_mult(c, a, b, KPRIMEC, LASTLIMBMASK, MPLIMB);
	}
	
	//gmp_printf("0x%Nx\n0x%Nx\n0x%Nx\n", a, GMPLIMB, b, GMPLIMB, c, GMPLIMB);
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randhlimbs(a, MPLIMB);
		randhlimbs(b, MPLIMB);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				gmp_mult(c, a, b, KPRIMEC, LASTLIMBMASK, MPLIMB);
				gmp_mult(b, c, a, KPRIMEC, LASTLIMBMASK, MPLIMB);
				gmp_mult(a, b, c, KPRIMEC, LASTLIMBMASK, MPLIMB);
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
	
	/*gmp_printf("0x%Nx\n0x%Nx\n0x%Nx\n", a, 9, b, 9, c, 9);
	gmp_printf("mtmp = 0x%Nx\n", mtmp, 18);
	gmp_printf("shift = 0x%Nx\n", shift, 9);*/
	
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
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				gmp_mult(c, a, b);
				gmp_mult(a, b, c);
				gmp_mult(b, c, a);
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

uint64_t do_pmersbench(void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]), const uint64_t W, const uint64_t NBLIMB, void (*langconvert)(uint64_t a[], uint64_t b[]))
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	uint64_t a[NBLIMB], b[NBLIMB], c[NBLIMB];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randhlimbs(a, NBLIMB);
		langconvert(a,a);
		randhlimbs(b, NBLIMB);
		langconvert(b,b);
		pmersmult(c, a, b);
	}
	
	/*bernhprint(a);
	bernhprint(b);
	bernhprint(c);*/
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randhlimbs(a, NBLIMB);
		langconvert(a,a);
		randhlimbs(b, NBLIMB);
		langconvert(b,b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint64_t soak=0; soak < W/3; soak++)
			{
				pmersmult(c, a, b);
				pmersmult(a, b, c);
				pmersmult(b, c, a);
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
		//medianTimer += timermin;
		free(statTimer);
	}
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

mp_limb_t tmp[20], carry, carry2[2];
void gmpmulmod2k(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB)
{
	mp_limb_t LASTLIMBMASK = (1ULL<<LASTLIMB)-1;
	mpn_mul_n(tmp, a, b, MPLIMB);
	carry = mpn_mul_1(c, tmp + MPLIMB, MPLIMB, KPRIMEC<<(64-LASTLIMB));
	carry += mpn_add_n(c, c, tmp, MPLIMB);
	carry = mpn_mul_1(carry2, &carry, 1, (KPRIMEC<<(64-LASTLIMB)));
	mpn_add_1(c, c, MPLIMB, carry2[0]);
	mpn_add_1(c + 1, c + 1, MPLIMB - 1, carry); 
	mpn_add_1(c, c, MPLIMB, ((c[MPLIMB-1] & (0xffffffffffffffff ^ LASTLIMBMASK))>>LASTLIMB) * KPRIMEC);
	c[MPLIMB-1] &= (LASTLIMBMASK);
}

void mersenne521(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	mpn_mul_n(mtmp, a, b, 9);//gmp_printf("\n0x%Nx\n", mtmp, 18);
	mpn_rshift(shift, mtmp + 8, 9, 9);//gmp_printf("0x%Nx\n\n", shift, 9);
	mtmp[8] &= 511;
	c[8] = 0;
	mpn_add_n(c, shift, mtmp, 9); //printf("c = 0x%lx\n", c[8]);
	c[0] += ((c[8] & 512)>0); c[8] &= 511;
}

void C25519convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[4];
	for(uint64_t i = 0; i < 4; i++)
		ax[i] = a[i];
	ax[4] &= 0x7fffffffffffffff;
	alang[0] = ax[0] % (1ULL<<51);
	alang[1] = ((ax[0] >> 51) | (ax[1] << 13)) % (1ULL<<51);
	alang[2] = ((ax[1] >> 38) | (ax[2] << 26)) % (1ULL<<51);
	alang[3] = ((ax[2] >> 25) | (ax[3] << 39)) % (1ULL<<51);
	alang[4] = ((ax[3] >> 12));
}

void M383convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[6];
	for(uint64_t i = 0; i < 6; i++)
		ax[i] = a[i];
	ax[5] &= 0x1fffffffffffffff;
	alang[0] = ax[0] % (1ULL<<55);
	alang[1] = ((ax[0] >> 55) | (ax[1] << 9)) % (1ULL<<55);
	alang[2] = ((ax[1] >> 46) | (ax[2] << 18)) % (1ULL<<55);
	alang[3] = ((ax[2] >> 37) | (ax[3] << 27)) % (1ULL<<55);
	alang[4] = ((ax[3] >> 28) | (ax[4] << 36)) % (1ULL<<55);
	alang[5] = ((ax[4] >> 19) | (ax[5] << 45)) % (1ULL<<55);
	alang[6] = ((ax[5] >> 10));
}

void C41417convert(uint64_t a[], uint64_t b[])
{
	uint64_t ax[7];
	for(uint64_t i = 0; i < 7; i++)
		ax[i] = a[i];
	ax[6] &= 0x3fffffff;
	b[0] = ax[0] % (1ULL<<60);
	b[1] = ((ax[0] >> 60) | (ax[1] << 4)) % (1ULL<<59);
	b[2] = ((ax[1] >> 55) | (ax[2] << 9)) % (1ULL<<59);
	b[3] = ((ax[2] >> 50) | (ax[3] << 14)) % (1ULL<<59);
	b[4] = ((ax[3] >> 45) | (ax[4] << 19)) % (1ULL<<59);
	b[5] = ((ax[4] >> 40) | (ax[5] << 24)) % (1ULL<<59);
	b[6] = ((ax[5] >> 35) | (ax[6] << 29)) % (1ULL<<59);
}

void Ed448convert(uint64_t a[], uint64_t b[])
{
	uint64_t ax[7];
	for(uint64_t i = 0; i < 7; i++)
		ax[i] = a[i];
	b[0] = ax[0] % (1ULL<<56);
	b[1] = ((ax[0] >> 56) | (ax[1] << 8)) % (1ULL<<56);
	b[2] = ((ax[1] >> 48) | (ax[2] << 16)) % (1ULL<<56);
	b[3] = ((ax[2] >> 40) | (ax[3] << 24)) % (1ULL<<56);
	b[4] = ((ax[3] >> 32) | (ax[4] << 32)) % (1ULL<<56);
	b[5] = ((ax[4] >> 24) | (ax[5] << 40)) % (1ULL<<56);
	b[6] = ((ax[5] >> 16) | (ax[6] << 48)) % (1ULL<<56);
	b[7] = ((ax[6] >> 8));
}

void M511convert(uint64_t a[], uint64_t b[])
{
	uint64_t ax[8];
	for(uint64_t i = 0; i < 8; i++)
		ax[i] = a[i];
	ax[7] &= 0x7fffffffffffffff;
	b[0] = ax[0] % (1ULL<<57);
	b[1] = ((ax[0] >> 57) | (ax[1] << 7)) % (1ULL<<57);
	b[2] = ((ax[1] >> 50) | (ax[2] << 14)) % (1ULL<<57);
	b[3] = ((ax[2] >> 43) | (ax[3] << 21)) % (1ULL<<57);
	b[4] = ((ax[3] >> 36) | (ax[4] << 28)) % (1ULL<<57);
	b[5] = ((ax[4] >> 29) | (ax[5] << 35)) % (1ULL<<57);
	b[6] = ((ax[5] >> 22) | (ax[6] << 42)) % (1ULL<<57);
	b[7] = ((ax[6] >> 15) | (ax[7] << 49)) % (1ULL<<57);
	b[8] = ((ax[7] >> 8));
}

void E521convert(uint64_t a[], uint64_t b[])
{
	uint64_t ax[9];
	for(uint64_t i = 0; i < 9; i++)
		ax[i] = a[i];
	ax[8] &= 0x1ff;
	b[0] = ax[0] % (1ULL<<58);
	b[1] = ((ax[0] >> 58) | (ax[1] << 6)) % (1ULL<<58);
	b[2] = ((ax[1] >> 52) | (ax[2] << 12)) % (1ULL<<58);
	b[3] = ((ax[2] >> 46) | (ax[3] << 18)) % (1ULL<<58);
	b[4] = ((ax[3] >> 40) | (ax[4] << 24)) % (1ULL<<58);
	b[5] = ((ax[4] >> 34) | (ax[5] << 30)) % (1ULL<<58);
	b[6] = ((ax[5] >> 28) | (ax[6] << 36)) % (1ULL<<58);
	b[7] = ((ax[6] >> 22) | (ax[7] << 42)) % (1ULL<<58);
	b[8] = ((ax[7] >> 16) | (ax[8] << 48));
}

uint64_t checkmodmul(const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB)
{
	mp_limb_t LASTLIMBMASK = (1ULL<<LASTLIMB)-1;
	uint64_t cpt = 0, tot = 0;
	mp_limb_t p[MPLIMB];
	for(int i = 0; i < MPLIMB; i++)
		p[i] = 0;
	p[MPLIMB - 1] = LASTLIMBMASK + 1;
	mpn_sub_1(p, p, MPLIMB, KPRIMEC);
	gmp_printf("0x%Nx\n", p, MPLIMB);
	mp_limb_t ax[MPLIMB],bx[MPLIMB],cx[MPLIMB],chk[MPLIMB], tmp[MPLIMB * 2], blank[MPLIMB*2];
	for(int i = 0; i < 1000000; i++)
	{
		printf("\b%ld\t%ld\r", cpt, cpt - (tot++));
		randhlimbs(ax, MPLIMB); ax[MPLIMB - 1] &= LASTLIMBMASK;
		randhlimbs(bx, MPLIMB); bx[MPLIMB - 1] &= LASTLIMBMASK;
		//gmp_printf("0x%Nx\n0x%Nx\n\n", ax, MPLIMB, bx, MPLIMB);
		mpn_mul_n(tmp, ax, bx, MPLIMB);
		mpn_tdiv_qr(blank, chk, 0, tmp, MPLIMB*2, p, MPLIMB);
		if(KPRIMEC == 1 && LASTLIMB == 9 && MPLIMB == 9)
			mersenne521(cx, ax, bx);
		else
			gmpmulmod2k(cx, ax, bx, KPRIMEC, LASTLIMB, MPLIMB);
		cpt += mpn_cmp(chk, cx, MPLIMB) == 0;
		//gmp_printf("0x%Nx\n0x%Nx\n", chk, MPLIMB, cx, MPLIMB);
	}
	return cpt;
}

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles;
	
	const uint64_t nbrepet = 303;
	
	
	
	
	cycles = do_bench(pmns_montg_mult, nbrepet);
	printf("pmns %ld\n", cycles);
	
	/*cycles = do_gmpbench(gmpmulmod2k, 19, 63, 4, nbrepet);
	printf("gmpmulmod2k 255-19 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 187, 63, 6, nbrepet);
	printf("gmpmulmod2k 383-187 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 17, 30, 7, nbrepet);
	printf("gmpmulmod2k 414-17 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 187, 63, 8, nbrepet);
	printf("gmpmulmod2k 511-187 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 1, 9, 9, nbrepet);
	printf("gmpmulmod2k 521-1 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 361, 63, 16, nbrepet);
	printf("gmpmulmod2k 1023-361 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 297, 18, 32, nbrepet);
	printf("gmpmulmod2k 2002-297 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 3819, 11, 63, nbrepet);
	printf("gmpmulmod2k 3979-3819 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 241, 5, 123, nbrepet);
	printf("gmpmulmod2k 7813-241 %ld\n", cycles);
	
	cycles = do_gmpbench(gmpmulmod2k, 2149, 1, 128, nbrepet);
	printf("gmpmulmod2k 8129-2149 %ld\n", cycles);*/
	
	/*printf("%ld\n", checkmodmul(19, 63, 4));
	printf("%ld\n", checkmodmul(187, 63, 6));
	printf("%ld\n", checkmodmul(17, 30, 7));
	printf("%ld\n", checkmodmul(187, 63, 8));*/
	//printf("%ld\n", checkmodmul(1, 9, 9));
	//printf("%ld\n", checkmodmul(361, 63, 16));
	//printf("%ld\n", checkmodmul(297, 18, 32));
	//printf("%ld\n", checkmodmul(3819, 11, 63));
	//printf("%ld\n", checkmodmul(241, 5, 123));
	//printf("%ld\n", checkmodmul(2149, 1, 128));
	
	
	/*cycles = do_mersennebench(mersenne521, nbrepet);
	printf("%ld\n", cycles);*/
	
	/*mp_limb_t plimbs[MPLIMB] = {0};
	plimbs[MPLIMB-1] = 0x8000000000000000;*/
	//mp_limb_t plimbs[MPLIMB] = {0};
	//plimbs[MPLIMB-1] = 0x40000000;
	//mpn_sub_1(plimbs, plimbs, MPLIMB, KPRIMEC);
	//mpn_sub_1(plimbs + 3, plimbs + 3, 4, 1ULL << 32);
	//mpn_sub_1(plimbs, plimbs, 8, 1);
	/*mp_limb_t p511[8] = {0};
	p511[7] = 0x8000000000000000;
	mpn_sub_1(p511, p511, 8, 187);
	mp_limb_t p521[9] = {0};
	p521[8] = 0x200;
	mpn_sub_1(p521, p521, 9, 1);
	
	
	//gmp_printf("0x%Nx\n", plimbs, MPLIMB);
	gmp_printf("0x%Nx\n", p511, 8);
	
	const uint64_t langlimbcnt = 9, mplimb = 8;
	uint64_t alang[langlimbcnt], blang[langlimbcnt], clang[langlimbcnt];
	mp_limb_t ax[mplimb],bx[mplimb],cx[mplimb],chk[mplimb], tmp[mplimb * 2], blank[mplimb*2], carry;
	uint64_t cpt = 0, tot=0;
	for(uint64_t i = 0; i < 1000; i++)
	{
		//randhlimbs(ax, MPLIMB);
		//randhlimbs(bx, MPLIMB);
		//ax[MPLIMB-1] &= 0x1fffffffffffffff;
		//bx[MPLIMB-1] &= 0x1fffffffffffffff;
		//ax[MPLIMB-1] &= 0x3fffffff;
		//bx[MPLIMB-1] &= 0x3fffffff;
		//mpn_mul_n(tmp, ax, bx, MPLIMB);
		//gmp_printf("0x%Nx\n", tmp, 2*MPLIMB);
		//mpn_tdiv_qr(blank, chk, 0, tmp, MPLIMB*2, plimbs, MPLIMB);
		//carry = mpn_mul_1(cx, tmp + MPLIMB, MPLIMB, KPRIMEC*2);
		//carry += mpn_add_n(cx, cx, tmp, MPLIMB);
		//cx[0] += (carry*2 + ((cx[MPLIMB-1] & 0x8000000000000000) > 0))*KPRIMEC;
		//cx[MPLIMB-1] &= 0x7fffffffffffffff;
		//gmp_printf("0x%Nx\n0x%Nx\n\n", cx, MPLIMB, chk, MPLIMB);
		//gmpmulmod2k(cx, ax, bx); cpt += mpn_cmp(chk, cx, MPLIMB) == 0;
		//alang[0] = ax[0] % (1ULL<<51);
		//alang[1] = ((ax[0] >> 51) + (ax[1] << 13)) % (1ULL<<51);
		//alang[2] = ((ax[1] >> 38) + (ax[2] << 26)) % (1ULL<<51);
		//alang[3] = ((ax[2] >> 25) + (ax[3] << 39)) % (1ULL<<51);
		//alang[4] = ((ax[3] >> 12));
		//gmp_printf("%Nx\n%lx %lx %lx %lx %lx\n", ax, 4, alang[4], alang[3], alang[2], alang[1], alang[0]);
		//blang[0] = bx[0] % (1ULL<<51);
		//blang[1] = ((bx[0] >> 51) + (bx[1] << 13)) % (1ULL<<51);
		//blang[2] = ((bx[1] >> 38) + (bx[2] << 26)) % (1ULL<<51);
		//blang[3] = ((bx[2] >> 25) + (bx[3] << 39)) % (1ULL<<51);
		//blang[4] = ((bx[3] >> 12));
		//multMod25519(clang, alang, blang);
		//printf("0x%lx\n0x%llx\n\n", clang[0], chk[0] % (1ULL<<51));
		//cpt += ((clang[0] & chk[0]) == clang[0]) && (((chk[3] >> 12) & clang[4]) == clang[4]);
		//alang[0] = ax[0] % (1ULL<<55);
		//alang[1] = ((ax[0] >> 55) + (ax[1] << 9)) % (1ULL<<55);
		//alang[2] = ((ax[1] >> 46) + (ax[2] << 18)) % (1ULL<<55);
		//alang[3] = ((ax[2] >> 37) + (ax[3] << 27)) % (1ULL<<55);
		//alang[4] = ((ax[3] >> 28) + (ax[4] << 36)) % (1ULL<<55);
		//alang[5] = ((ax[4] >> 19) + (ax[5] << 45)) % (1ULL<<55);
		//alang[6] = ((ax[5] >> 10));
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx\n", ax, 6, alang[6], alang[5], alang[4], alang[3], alang[2], alang[1], alang[0]);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx\n", bx, 6, alang[6], alang[5], alang[4], alang[3], alang[2], alang[1], alang[0]);
		//blang[0] = bx[0] % (1ULL<<55);
		//blang[1] = ((bx[0] >> 55) + (bx[1] << 9)) % (1ULL<<55);
		//blang[2] = ((bx[1] >> 46) + (bx[2] << 18)) % (1ULL<<55);
		//blang[3] = ((bx[2] >> 37) + (bx[3] << 27)) % (1ULL<<55);
		//blang[4] = ((bx[3] >> 28) + (bx[4] << 36)) % (1ULL<<55);
		//blang[5] = ((bx[4] >> 19) + (bx[5] << 45)) % (1ULL<<55);
		//blang[6] = ((bx[5] >> 10));
		//multModM383(clang, alang, blang);
		//gmp_printf("0x%Nx\n0x%Nx\n0x%Nx\n", ax, MPLIMB, bx, MPLIMB, chk, MPLIMB);
		//printf("0x%lx\n0x%lx\n0x%llx\n\n", clang[6], clang[0], chk[0] % (1ULL<<55));
		//cpt += (((clang[0] & chk[0]) == clang[0]) &&
		//       ((((chk[0]>>55) | (chk[1]<<9)) & clang[1]) == clang[1]) &&
		//       ((((chk[1]>>46) | (chk[2]<<18)) & clang[2]) == clang[2]) &&
		//       ((((chk[2]>>37) | (chk[3]<<27)) & clang[3]) == clang[3]) &&
		//       ((((chk[3]>>28) | (chk[4]<<36)) & clang[4]) == clang[4]) &&
		//       ((((chk[4]>>19) | (chk[5]<<45)) & clang[5]) == clang[5]) &&
		//       ((((chk[5] >> 10) & clang[6]) == clang[6])));
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx\n", chk, 6, clang[6], clang[5], clang[4], clang[3], clang[2], clang[1], clang[0]);
		//alang[0] = ax[0] % (1ULL<<60);
		//alang[1] = ((ax[0] >> 60) + (ax[1] << 4)) % (1ULL<<59);
		//alang[2] = ((ax[1] >> 55) + (ax[2] << 9)) % (1ULL<<59);
		//alang[3] = ((ax[2] >> 50) + (ax[3] << 14)) % (1ULL<<59);
		//alang[4] = ((ax[3] >> 45) + (ax[4] << 19)) % (1ULL<<59);
		//alang[5] = ((ax[4] >> 40) + (ax[5] << 24)) % (1ULL<<59);
		//alang[6] = ((ax[5] >> 35) + (ax[6] << 29)) % (1ULL<<59);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx\n", ax, 7, alang[6], alang[5], alang[4], alang[3], alang[2], alang[1], alang[0]);
		//blang[0] = bx[0] % (1ULL<<60);
		//blang[1] = ((bx[0] >> 60) + (bx[1] << 4)) % (1ULL<<59);
		//blang[2] = ((bx[1] >> 55) + (bx[2] << 9)) % (1ULL<<59);
		//blang[3] = ((bx[2] >> 50) + (bx[3] << 14)) % (1ULL<<59);
		//blang[4] = ((bx[3] >> 45) + (bx[4] << 19)) % (1ULL<<59);
		//blang[5] = ((bx[4] >> 40) + (bx[5] << 24)) % (1ULL<<59);
		//blang[6] = ((bx[5] >> 35) + (bx[6] << 29)) % (1ULL<<59);
		//multModC41417(clang, alang, blang);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx\n", chk, 7, clang[6], clang[5], clang[4], clang[3], clang[2], clang[1], clang[0]);
		//cpt += (((clang[0] & chk[0]) == clang[0]) &&
		//       ((((chk[0]>>60) | (chk[1]<<4)) & clang[1]) == clang[1]) &&
		//       ((((chk[1]>>55) | (chk[2]<<9)) & clang[2]) == clang[2]) &&
		//       ((((chk[2]>>50) | (chk[3]<<14)) & clang[3]) == clang[3]) &&
		//       ((((chk[3]>>45) | (chk[4]<<19)) & clang[4]) == clang[4]) &&
		//       ((((chk[4]>>40) | (chk[5]<<24)) & clang[5]) == clang[5]) &&
		//       ((((chk[5]>>35) | (chk[6]<<29)) & clang[6]) == clang[6]));
		//printf("\b%ld\t%ld\r", cpt, cpt - (tot++));
		//Ed448convert(ax, alang);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx %lx\n", ax, 7, alang[7], alang[6], alang[5], alang[4], alang[3], alang[2], alang[1], alang[0]);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx %lx\n", bx, 7, alang[7], alang[6], alang[5], alang[4], alang[3], alang[2], alang[1], alang[0]);
		//Ed448convert(bx, blang);
		//multModEd448(clang, alang, blang);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx %lx\n", chk, 7, clang[7], clang[6], clang[5], clang[4], clang[3], clang[2], clang[1], clang[0]);
		//cpt += (((clang[0] & chk[0]) == clang[0]) &&
		//       ((((chk[0]>>56) | (chk[1]<<8)) & clang[1]) == clang[1]) &&
		//       ((((chk[1]>>48) | (chk[2]<<16)) & clang[2]) == clang[2]) &&
		//       ((((chk[2]>>40) | (chk[3]<<24)) & clang[3]) == clang[3]) &&
		//       ((((chk[3]>>32) | (chk[4]<<32)) & clang[4]) == clang[4]) &&
		//       ((((chk[4]>>24) | (chk[5]<<40)) & clang[5]) == clang[5]) &&
		//       ((((chk[5]>>16) | (chk[6]<<48)) & clang[6]) == clang[6]) &&
		//       (((chk[6] >> 8) & clang[7]) == clang[7]));
		//randhlimbs(ax, 8); ax[7] &= 0x7fffffffffffffff;
		//randhlimbs(bx, 8); bx[7] &= 0x7fffffffffffffff;
		//mpn_mul_n(tmp, ax, bx, 8);
		//mpn_tdiv_qr(blank, chk, 0, tmp, 8*2, p511, 8);
		//M511convert(ax, alang);
		//M511convert(bx, blang);
		//multModM511(clang, alang, blang);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx %lx %lx\n", chk, 8, clang[8], clang[7], clang[6], clang[5], clang[4], clang[3], clang[2], clang[1], clang[0]);
		//cpt += (((clang[0] & chk[0]) == clang[0]) &&
		//       ((((chk[0]>>57) | (chk[1]<<7)) & clang[1]) == clang[1]) &&
		//       ((((chk[1]>>50) | (chk[2]<<14)) & clang[2]) == clang[2]) &&
		//       ((((chk[2]>>43) | (chk[3]<<21)) & clang[3]) == clang[3]) &&
		//       ((((chk[3]>>36) | (chk[4]<<28)) & clang[4]) == clang[4]) &&
		//       ((((chk[4]>>29) | (chk[5]<<35)) & clang[5]) == clang[5]) &&
		//       ((((chk[5]>>22) | (chk[6]<<42)) & clang[6]) == clang[6]) &&
		//       ((((chk[6]>>15) | (chk[7]<<49)) & clang[7]) == clang[7]) &&
		//       (((chk[7] >> 8) & clang[8]) == clang[8]));
		//randhlimbs(ax, 9); ax[8] &= 0x1ff;
		//randhlimbs(bx, 9); bx[8] &= 0x1ff;
		//mpn_mul_n(tmp, ax, bx, 9);
		//mpn_tdiv_qr(blank, chk, 0, tmp, 9*2, p521, 9);
		//E521convert(ax, alang);
		//E521convert(bx, blang);
		//multModE521(clang, alang, blang);
		//gmp_printf("0x%Nx\n%lx %lx %lx %lx %lx %lx %lx %lx %lx\n", chk, 9, clang[8], clang[7], clang[6], clang[5], clang[4], clang[3], clang[2], clang[1], clang[0]);
		//cpt += (((clang[0] & chk[0]) == clang[0]) &&
		//       ((((chk[0]>>58) | (chk[1]<<6)) & clang[1]) == clang[1]) &&
		//       ((((chk[1]>>52) | (chk[2]<<12)) & clang[2]) == clang[2]) &&
		//       ((((chk[2]>>46) | (chk[3]<<18)) & clang[3]) == clang[3]) &&
		//       ((((chk[3]>>40) | (chk[4]<<24)) & clang[4]) == clang[4]) &&
		//       ((((chk[4]>>34) | (chk[5]<<30)) & clang[5]) == clang[5]) &&
		//       ((((chk[5]>>28) | (chk[6]<<36)) & clang[6]) == clang[6]) &&
		//       ((((chk[6]>>22) | (chk[7]<<42)) & clang[7]) == clang[7]) &&
		//       ((((chk[7]>>16) | (chk[8]<<48)) & clang[8]) == clang[8]));
	}
	printf("%ld\n", cpt);*/
	
	
	
	/*cycles = do_pmersbench(multMod25519, nbrepet, 5, C25519convert);
	printf("C25519 %ld\n", cycles);
	
	cycles = do_pmersbench(multModM383, nbrepet, 7, M383convert);
	printf("M-383 %ld\n", cycles);
	
	cycles = do_pmersbench(multModC41417, nbrepet, 7, C41417convert);
	printf("Curve41417 %ld\n", cycles);
	
	cycles = do_pmersbench(multModEd448, nbrepet, 8, Ed448convert);
	printf("Ed448-Goldilocks %ld\n", cycles);
	
	cycles = do_pmersbench(multModM511, nbrepet, 9, M511convert);
	printf("M-511 %ld\n", cycles);
	
	cycles = do_pmersbench(multModE521, nbrepet, 9, E521convert);
	printf("E-521 %ld\n", cycles);*/
	
	return 0;
}


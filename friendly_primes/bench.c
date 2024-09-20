#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "bench.h"
#include "eccoptimizedcode.h"

#define NTEST 501
#define NSAMPLES 1001

extern uint8_t ALPHA;

int64_t randomint64(void)
{
	// Function that generates a random integer on 64 bits.
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}


void randpoly(poly P, uint64_t RHO)
{
	// Function that generates a random polynomial with all coefficients
	// absolutely bounded by RHO.
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = (randomint64() % RHO) * (1 + (rand() & 1) * -2);
}

void randhlimbs(mp_limb_t g[], const uint64_t limbs)
{
	// Function that generates a random limb array.
	for(uint64_t i = 0; i < limbs; i++)
		g[i] = randomint64();
}

void randmersennelimbs(mp_limb_t g[])
{
	// Function that generates a random number for E-521.
	randhlimbs(g, 9);
	g[8] = g[8] % 512;
}


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

/**** From here on, wrapper functions.
*****/

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint16_t N, const uint64_t RHO, const uint64_t W)
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
			printf("\b%d\t%d\r", i, j);
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
	
	printf("                                          \r");
	
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
			printf("\b%d\t%d\r", i, j);
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
	
	printf("                                          \r");
	
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
			printf("\b%d\t%d\r", i, j);
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
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

uint64_t do_pmersbench(void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]), const uint64_t W, const uint64_t NBLIMB, void (*convertfunc)(uint64_t a[], uint64_t b[]))
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	uint64_t a[NBLIMB], b[NBLIMB], c[NBLIMB];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randhlimbs(a, NBLIMB);
		convertfunc(a,a);
		randhlimbs(b, NBLIMB);
		convertfunc(b,b);
		pmersmult(c, a, b);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randhlimbs(a, NBLIMB);
		convertfunc(a,a);
		randhlimbs(b, NBLIMB);
		convertfunc(b,b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			printf("\b%d\t%d\r", i, j);
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
		free(statTimer);
	}
	
	printf("                                          \r");
	
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

mp_limb_t tmp[70], carry, carry2[2];
void gmpmulmod2k(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB)
{
	// Adapted GMP modular multiplication for pseudo mersenne modulus.
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

mp_limb_t tmpbig[270];
void gmpmulmod2kbig(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB)
{
	// Adapted GMP modular multiplication for big pseudo mersenne modulus.
	mp_limb_t LASTLIMBMASK = (1ULL<<LASTLIMB)-1;
	mpn_mul_n(tmpbig, a, b, MPLIMB);
	carry = mpn_mul_1(c, tmpbig + MPLIMB, MPLIMB, KPRIMEC)<<(64-LASTLIMB);
	carry += mpn_lshift(c, c,  MPLIMB, 64-LASTLIMB);
	carry += mpn_add_n(c, c, tmpbig, MPLIMB);
	carry = mpn_mul_1(carry2, &carry, 1, KPRIMEC)<<(64-LASTLIMB);
	carry += mpn_lshift(carry2, carry2, 2, 64-LASTLIMB);
	mpn_add_1(c, c, MPLIMB, carry2[0]);
	carry += carry2[1];
	mpn_add_1(c + 1, c + 1, MPLIMB - 1, carry);
	unsigned __int128 carry3 = (__int128)((c[MPLIMB-1] & (0xffffffffffffffff ^ LASTLIMBMASK))>>LASTLIMB) * KPRIMEC;
	mpn_add_1(c, c, MPLIMB, (uint64_t) carry3);
	mpn_add_1(c + 1, c + 1, MPLIMB - 1, (uint64_t) (carry3>>64));
	c[MPLIMB-1] &= (LASTLIMBMASK);
}

void mersenne521(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[])
{
	// Adapted GMP modular multiplication for 2^521 - 1.
	mpn_mul_n(mtmp, a, b, 9);//gmp_printf("\n0x%Nx\n", mtmp, 18);
	mpn_rshift(shift, mtmp + 8, 9, 9);//gmp_printf("0x%Nx\n\n", shift, 9);
	mtmp[8] &= 511;
	c[8] = 0;
	mpn_add_n(c, shift, mtmp, 9); //printf("c = 0x%lx\n", c[8]);
	c[0] += ((c[8] & 512)>0); c[8] &= 511;
}

/**** From here on, conversion functions.
****/

void PoC244convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[4];
	for(uint64_t i = 0; i < 4; i++)
		ax[i] = a[i];
	ax[3] &= 0xfffffffffffff;
	alang[0] = ax[0] % (1ULL<<49);
	alang[1] = ((ax[0] >> 49) | (ax[1] << 15)) % (1ULL<<49);
	alang[2] = ((ax[1] >> 34) | (ax[2] << 30)) % (1ULL<<49);
	alang[3] = ((ax[2] >> 19) | (ax[3] << 45)) % (1ULL<<49);
	alang[4] = ((ax[3] >> 4));
}

void PoC244convertv2(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[4];
	for(uint64_t i = 0; i < 4; i++)
		ax[i] = a[i];
	ax[3] &= 0xfffffffffffff;
	alang[0] = ax[0] % (1ULL<<61);
	alang[1] = ((ax[0] >> 61) | (ax[1] << 3)) % (1ULL<<61);
	alang[2] = ((ax[1] >> 58) | (ax[2] << 6)) % (1ULL<<61);
	alang[3] = ((ax[2] >> 55) | (ax[3] << 9)) % (1ULL<<61);
}

void PoC297convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[5];
	for(uint64_t i = 0; i < 5; i++)
		ax[i] = a[i];
	ax[4] &= 0x1ffffffffff;
	alang[0] = ax[0] % (1ULL<<50);
	alang[1] = ((ax[0] >> 50) | (ax[1] << 14)) % (1ULL<<50);
	alang[2] = ((ax[1] >> 36) | (ax[2] << 28)) % (1ULL<<50);
	alang[3] = ((ax[2] >> 22) | (ax[3] << 42)) % (1ULL<<50);
	alang[4] = ((ax[3] >> 8) | (ax[4] << 56)) % (1ULL<<50);
	alang[5] = ((ax[3] >> 58) | (ax[4] << 6)) % (1ULL<<47);
}

void PoC354convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[6];
	for(uint64_t i = 0; i < 6; i++)
		ax[i] = a[i];
	ax[5] &= 0x3ffffffff;
	alang[0] = ax[0] % (1ULL<<51);
	alang[1] = ((ax[0] >> 51) | (ax[1] << 13)) % (1ULL<<51);
	alang[2] = ((ax[1] >> 38) | (ax[2] << 26)) % (1ULL<<51);
	alang[3] = ((ax[2] >> 25) | (ax[3] << 39)) % (1ULL<<51);
	alang[4] = ((ax[3] >> 12) | (ax[4] << 52)) % (1ULL<<50);
	alang[5] = ((ax[3] >> 62) | (ax[4] << 2)) % (1ULL<<50);
	alang[6] = ((ax[4] >> 48) | (ax[5] << 16)) % (1ULL<<50);
}

void PoC480convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[8];
	for(uint64_t i = 0; i < 8; i++)
		ax[i] = a[i];
	ax[7] &= 0xffffffff;
	alang[0] = ax[0] % (1ULL<<54);
	alang[1] = ((ax[0] >> 54) | (ax[1] << 10)) % (1ULL<<54);
	alang[2] = ((ax[1] >> 44) | (ax[2] << 20)) % (1ULL<<54);
	alang[3] = ((ax[2] >> 34) | (ax[3] << 30)) % (1ULL<<53);
	alang[4] = ((ax[3] >> 23) | (ax[4] << 41)) % (1ULL<<53);
	alang[5] = ((ax[4] >> 12) | (ax[5] << 52)) % (1ULL<<53);
	alang[6] = ((ax[5] >> 1) | (ax[6] << 63)) % (1ULL<<53);
	alang[7] = ((ax[5] >> 54) | (ax[6] << 10)) % (1ULL<<53);
	alang[8] = ((ax[6] >> 43) | (ax[7] << 21)) % (1ULL<<53);
	alang[9] = ((ax[7] >> 32) % (1ULL<<53));
}

void C25519convert(uint64_t a[], uint64_t alang[])
{
	uint64_t ax[4];
	for(uint64_t i = 0; i < 4; i++)
		ax[i] = a[i];
	ax[3] &= 0x7fffffffffffffff;
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
	ax[5] &= 0x7fffffffffffffff;
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
	// Function that checks 1000000 random multiplications for errors.
	mp_limb_t LASTLIMBMASK = (1ULL<<LASTLIMB)-1;
	uint64_t cpt = 0;
	mp_limb_t p[MPLIMB];
	for(uint i = 0; i < MPLIMB; i++)
		p[i] = 0;
	p[MPLIMB - 1] = LASTLIMBMASK + 1;
	mpn_sub_1(p, p, MPLIMB, KPRIMEC);
	mp_limb_t ax[MPLIMB],bx[MPLIMB],cx[MPLIMB],chk[MPLIMB], tmp[MPLIMB * 2], blank[MPLIMB*2];
	
	uint64_t NBLOOPS = LOOPCHKNMB;
	if(MPLIMB >= 63)
		NBLOOPS /= 10;  // otherwise it takes a while
	
	for(uint64_t i = 0; i < NBLOOPS; i++)
	{
		printf("\b%ld\t%ld\r", i, cpt);
		randhlimbs(ax, MPLIMB); ax[MPLIMB - 1] &= LASTLIMBMASK;
		randhlimbs(bx, MPLIMB); bx[MPLIMB - 1] &= LASTLIMBMASK;
		mpn_mul_n(tmp, ax, bx, MPLIMB);
		mpn_tdiv_qr(blank, chk, 0, tmp, MPLIMB*2, p, MPLIMB);
		if(KPRIMEC == 1 && LASTLIMB == 9 && MPLIMB == 9)
			mersenne521(cx, ax, bx);
		else if(MPLIMB<63)
			gmpmulmod2k(cx, ax, bx, KPRIMEC, LASTLIMB, MPLIMB);
		else
			gmpmulmod2kbig(cx, ax, bx, KPRIMEC, LASTLIMB, MPLIMB);
		cpt += mpn_cmp(chk, cx, MPLIMB) == 0;
	}
	printf("                                          \r");
	return cpt;
}

int64_t _checkOpti(const uint16_t SIZEP, mp_limb_t pmers[], const uint16_t NSIZE, void (*convertfunc)(uint64_t a[], uint64_t b[]), void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]))
{
	// Function that checks 1000000 random multiplications for errors.
	uint64_t alang[NSIZE], blang[NSIZE], clang[NSIZE], chklang[NSIZE];
	mp_limb_t ax[SIZEP],bx[SIZEP],chk[SIZEP], tmp[SIZEP*2], blank[SIZEP*2];
	ax[SIZEP - 1] = 0; // If I don't do this gcc thinks it's uninitialized.
	int64_t cpt = 0;
	_Bool check;
	
	for(int i = 0; i < LOOPCHKNMB; i++)
	{
		printf("\b%d\t%ld\r", i, cpt);
		randhlimbs(tmp, SIZEP);
		mpn_tdiv_qr(blank, ax, 0, tmp, SIZEP, pmers, SIZEP);
		randhlimbs(tmp, SIZEP);
		mpn_tdiv_qr(blank, bx, 0, tmp, SIZEP, pmers, SIZEP);
		mpn_mul_n(tmp, ax, bx, SIZEP);
		mpn_tdiv_qr(blank, chk, 0, tmp, SIZEP*2, pmers, SIZEP);
		convertfunc(ax, alang);
		convertfunc(bx, blang);
		pmersmult(clang, alang, blang);
		convertfunc(chk, chklang);
		check = 1;
		for(int j = 0; j < NSIZE; j++)
			check = check && (chklang[j] == clang[j]);
		cpt += check;
	}
	printf("                                          \r");
	return cpt;
}

static inline int64_t checkGoldilocks(void (*multModEd448)(uint64_t c[], uint64_t a[], uint64_t b[]))
{
	// Wrapper function for Goldilocks.
	mp_limb_t p448[8];
	for(int i = 0; i < 7; i++)
		p448[i] = 0;
	p448[7] = 0x1;
	mpn_sub_1(p448+3, p448+3, 5, 0x100000000);
	mpn_sub_1(p448,p448,7,1);
	return _checkOpti(7, p448, 8, Ed448convert, multModEd448);
}

int64_t checkOpti(const uint16_t SIZEP, const uint64_t PMASK, const uint64_t PSUB, const uint16_t NSIZE, void (*convertfunc)(uint64_t a[], uint64_t b[]), void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]))
{
	// Wrapper function for Mersenne and Pseudo Mersenne.
	
	if(SIZEP == 0 && PMASK == 0 && PSUB == 0 && NSIZE == 0 && convertfunc == 0)
		return checkGoldilocks(pmersmult);
	else
	{
		mp_limb_t pmers[SIZEP];
		for(int i = 0; i < SIZEP - 1; i++)
			pmers[i] = 0;
		pmers[SIZEP - 1] = PMASK;
		mpn_sub_1(pmers, pmers, SIZEP, PSUB);
		return _checkOpti(SIZEP, pmers, NSIZE, convertfunc, pmersmult);
	}
}


void horner_eval(mpz_t res, const poly P, const int64_t gamma)
{
	// Function that evaluates the polynomial P in gamma.
	// Uses the horner algorithm.
	
	mpz_t A;
	mpz_init(A);
	
	uint16_t N = P->deg;
	mpz_set_si(A, P->t[N-1]);
	
	for(int i = N - 2; i >= 0; i--)
	{
		mpz_mul_ui(A, A, gamma);
		if(P->t[i] > 0)
			mpz_add_ui(A, A, P->t[i]);
		else
			mpz_sub_ui(A, A, -P->t[i]);
	}
	
	mpz_set(res, A);
	
	mpz_clear(A);
}

int64_t checkPmns(const uint16_t N, const uint64_t RHO, const uint64_t GAMMA, mp_limb_t prime[], void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly))
{
	_poly a, b, c;
	int64_t atab[N], btab[N], ctab[N], cpt = 0;
	a.deg = N; b.deg = N; c.deg = N;
	a.t = atab; b.t = btab; c.t = ctab;
	mpz_t A, B, C, CC, PHI, PRIME;
	mpz_inits(A, B, C, CC, PHI, PRIME, NULL);
	mp_limb_t phitab[2] = {0, 1};
	PHI->_mp_size = 2;
	PHI->_mp_alloc = 2;
	PHI->_mp_d = phitab;
	PRIME->_mp_size = N;
	PRIME->_mp_alloc = N;
	PRIME->_mp_d = prime;
	while(PRIME->_mp_d[PRIME->_mp_size - 1] == 0)
		PRIME->_mp_size--;
	
	uint64_t morethanrho = 0;
	uint64_t NBLOOPS = LOOPCHKNMB;
	if(N >= 72)
		NBLOOPS /= 10;  // otherwise it takes a while
	for(uint64_t i = 0; i < NBLOOPS; i++)
	{
		printf("\b%ld\t%ld\r", i, cpt);
		randpoly(&a, RHO);
		randpoly(&b, RHO);
		pmns_mult(&c, &a, &b);
		
		for(int j = 0; j < N; j++)
			morethanrho += ((c.t[j] > (int64_t) RHO) || (c.t[j] < -((int64_t) RHO)));
		
		horner_eval(A, &a, GAMMA);
		horner_eval(B, &b, GAMMA);
		horner_eval(C, &c, GAMMA);
		mpz_mul(C, C, PHI);
		
		mpz_mod(A, A, PRIME);
		mpz_mod(B, B, PRIME);
		mpz_mod(C, C, PRIME);
		
		mpz_mul(CC, A, B);
		mpz_mul_ui(CC, CC, ALPHA);
		mpz_mod(CC, CC, PRIME);
		
		cpt += mpz_cmp(C, CC) == 0;
	}
	printf("                                          \r");
	
	if(morethanrho)
		printf("More than RHO: %lu\n", morethanrho);
	
	PRIME->_mp_size = 0;
	PRIME->_mp_alloc = 0;
	PRIME->_mp_d = NULL;
	PHI->_mp_size = 0;
	PHI->_mp_alloc = 0;
	PHI->_mp_d = NULL;
	mpz_clears(A, B, C, CC, PHI, PRIME, NULL);
	
	return cpt;
}

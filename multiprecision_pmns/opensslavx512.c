#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#include <gmp.h>
#include <openssl/bn.h>
#include <assert.h>

#define NTEST 501
#define NSAMPLES 1001
#define NBREPET 33

uint32_t OPENSSL_ia32cap_P[8] = { 0 };

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

extern void ossl_rsaz_amm52x20_x1_ifma256(BN_ULONG *res,
                                    const BN_ULONG *a,
                                    const BN_ULONG *b,
                                    const BN_ULONG *m,
                                    BN_ULONG k0);

extern void ossl_rsaz_amm52x30_x1_ifma256(BN_ULONG *res,
                                    const BN_ULONG *a,
                                    const BN_ULONG *b,
                                    const BN_ULONG *m,
                                    BN_ULONG k0);

extern void ossl_rsaz_amm52x40_x1_ifma256(BN_ULONG *res,
                                    const BN_ULONG *a,
                                    const BN_ULONG *b,
                                    const BN_ULONG *m,
                                    BN_ULONG k0);

typedef struct bignum_st BIGNUM;

struct bignum_st
       {
       BN_ULONG *d;    /* Pointer to an array of 'BN_BITS2' bit chunks. */
       int top;        /* Index of last used d +1. */
       /* The next are internal book keeping for bn_expand. */
       int dmax;       /* Size of the d array. */
       int neg;        /* one if the number is negative */
       int flags;
       };

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

/****

	End of section.

*****/

uint64_t opensslbenchavx512(BN_ULONG* P, BN_ULONG invP, uint16_t W, void (*openssl_mul)(BN_ULONG *res, const BN_ULONG *a, const BN_ULONG *b, const BN_ULONG *m, BN_ULONG k0))
{
	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1, t2, diff_t, *statTimer;
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	BN_ULONG a[LOG2P/52 + 1], b[LOG2P/52 + 1], c[LOG2P/52 + 1];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		for(int k = 0; k < LOG2P/52; k++)
			a[k] = randomint64() % (1ULL << 52);
		for(int k = 0; k < LOG2P/52; k++)
			b[k] = randomint64() % (1ULL << 52);
		a[LOG2P/52] = randomint64() % P[LOG2P/52];
		b[LOG2P/52] = randomint64() % P[LOG2P/52];
		openssl_mul(c, a, b, P, invP);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		for(int k = 0; k < LOG2P/52; k++)
			a[k] = randomint64() % (1ULL << 52);
		for(int k = 0; k < LOG2P/52; k++)
			b[k] = randomint64() % (1ULL << 52);
		a[LOG2P/52] = randomint64() % P[LOG2P/52];
		b[LOG2P/52] = randomint64() % P[LOG2P/52];		
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W/3; soak++)
			{
				openssl_mul(c, a, b, P, invP);
				openssl_mul(a, b, c, P, invP);
				openssl_mul(b, c, a, P, invP);
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

int main(void)
{
	mpz_t A, B, P, R, invP;
	mpz_inits(A, B, P, R, invP, NULL);
	
	mpz_set_ui(R, 1);
	
	mpz_mul_2exp(R, R, LOG2P);
	
	long seed = time(NULL);
	srand((unsigned) (time(&seed)));
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	mpz_urandomb(P, r, LOG2P);
	
	P->_mp_d[LOG2P/64-1] |= 1ULL << 63;
	
	mpz_nextprime(P, P);
	
	mpz_invert(invP, P, R);
	
	/*gmp_printf("\n   P =  0x%Zx\n", P);
	gmp_printf("invP =  0x%Zx\n", invP);
	gmp_printf("   R = 0x%Zx\n", R);*/
	
	BN_CTX *ctx = BN_CTX_new();
	BN_CTX_start(ctx);
	
	BN_ULONG p[2*LOG2P/52 + 1] = {0}, invp;
	
	for(int i = 0; i < LOG2P/52 + 1; i++)
	{
		p[i] = P->_mp_d[0] % (1ULL<<52);
		mpz_tdiv_r_2exp(P, P, 52);
	}
	
	invp = (-invP->_mp_d[0]) % (1ULL<<52);
	
	#if LOG2P == 1024
	printf("OpenSSL Montgomery CIOS: %ld\n", opensslbenchavx512(p, invp, NBREPET, ossl_rsaz_amm52x20_x1_ifma256));
	#endif
	#if LOG2P == 1536
	printf("OpenSSL Montgomery CIOS: %ld\n", opensslbenchavx512(p, invp, NBREPET, ossl_rsaz_amm52x30_x1_ifma256));
	#endif
	#if LOG2P == 2048
	printf("OpenSSL Montgomery CIOS: %ld\n", opensslbenchavx512(p, invp, NBREPET, ossl_rsaz_amm52x40_x1_ifma256));
	#endif
	
	mpz_clears(A, B, P, R, invP, NULL);
	gmp_randclear(r);
	
	BN_CTX_end(ctx);
	BN_CTX_free(ctx);
	
	return 0;
}

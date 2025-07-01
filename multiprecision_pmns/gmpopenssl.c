#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
#include <openssl/bn.h>

#define UNUSED(X) (void)(X)

#define NTEST 501
#define NSAMPLES 1001

#define mpn_redc_1 __MPN(redc_1)
__GMP_DECLSPEC mp_limb_t mpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);

#define mpn_redc_n __MPN(redc_n)
__GMP_DECLSPEC void mpn_redc_n (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);

#define	 mpn_binvert __MPN(binvert)
__GMP_DECLSPEC void mpn_binvert (mp_ptr, mp_srcptr, mp_size_t, mp_ptr);

#define mpn_mullo_n __MPN(mullo_n)
__GMP_DECLSPEC void mpn_mullo_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);


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


void mpn_mont_mul_red_1(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2, const mp_limb_t *p_limbs, mp_limb_t mip0, int nb_limbs)
{
	mp_limb_t tmp_limbs[2*nb_limbs];
	
	mpn_mul_n(tmp_limbs, op1, op2, nb_limbs);
	
	if (mpn_redc_1(rop, tmp_limbs, p_limbs, nb_limbs, mip0) != 0)
		mpn_sub_n(rop, rop, p_limbs, nb_limbs);
}

void mpn_mont_mul_red_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2, const mp_limb_t *p_limbs, mp_limb_t *mip_limbs, int nb_limbs)
{
	mp_limb_t tmp_limbs[2*nb_limbs];
	
	mpn_mul_n(tmp_limbs, op1, op2, nb_limbs);
	mpn_redc_n(rop, tmp_limbs, p_limbs, nb_limbs, mip_limbs);
}

void gmp_montgomery_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *mip_limbs, mp_limb_t *soak, int nb_limbs)
{
	UNUSED(soak);
	mpn_mont_mul_red_n(a_limbs, a_limbs, b_limbs, p_limbs, mip_limbs, nb_limbs);
}

void gmp_montgomeryCIOS_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *mip0, mp_limb_t *soak, int nb_limbs)
{
	UNUSED(soak);
	mpn_mont_mul_red_1(a_limbs, a_limbs, b_limbs, p_limbs, *mip0, nb_limbs);
}

uint64_t gmpbench(mpz_t A, mpz_t B, mpz_t modul_p, gmp_randstate_t r, mp_limb_t *c_limbs, mp_limb_t *q_limbs, uint16_t W, void (*gmp_wrapper)(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
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

uint64_t opensslbench(mpz_t A, mpz_t B, mpz_t P, gmp_randstate_t r, BN_CTX* ctx, BN_MONT_CTX* mont_ctx, uint16_t W, int (*openssl_mul)(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, BN_MONT_CTX *mont, BN_CTX *ctx))
{
	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1, t2, diff_t, *statTimer;
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	char* stmp;
	
	
	BIGNUM* a = BN_CTX_get(ctx);
	BIGNUM* b = BN_CTX_get(ctx);
	BIGNUM* c = BN_CTX_get(ctx);
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		mpz_urandomm(A, r, P);
		mpz_urandomm(B, r, P);
		stmp = mpz_get_str(NULL, 10, A);
		BN_dec2bn(&a, stmp);
		free(stmp);
		stmp = mpz_get_str(NULL, 10, B);
		BN_dec2bn(&b, stmp);
		free(stmp);
		openssl_mul(c, a, b, mont_ctx, ctx);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		mpz_urandomm(A, r, P);
		mpz_urandomm(B, r, P);
		
		stmp = mpz_get_str(NULL, 10, A);
		BN_dec2bn(&a, stmp);
		free(stmp);
		stmp = mpz_get_str(NULL, 10, B);
		BN_dec2bn(&b, stmp);
		free(stmp);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W/3; soak++)
			{
				openssl_mul(c, a, b, mont_ctx, ctx);
				openssl_mul(a, b, c, mont_ctx, ctx);
				openssl_mul(b, c, a, mont_ctx, ctx);
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
	#if LOG2P >= 16*64
	const uint64_t W = 3;
	#else
	const uint64_t W = 33;
	#endif
	mpz_t P, R, invP;
	mpz_inits(P, R, invP, NULL);
	
	mpz_set_ui(R, 1);
	
	mpz_mul_2exp(R, R, LOG2P);
	
	unsigned long seed = time(NULL);
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
	
	mpz_t A, B, C;
	mpz_inits(A, B, C, NULL);
	
	//mp_limb_t c_limbs[LOG2P/32] = {0};
	mp_limb_t q_limbs[LOG2P/64+1] = {0};
	
	mp_limb_t mip0;
	mip0 = invP->_mp_d[0];
	mip0 = -mip0;
	
	mpz_urandomm(A, r, P);
	mpz_urandomm(B, r, P);
	mpz_set(C, A);
	
	/*gmp_printf("\nA = 0x%Zx\n", A);
	gmp_printf("B = 0x%Zx\n", B);
	gmp_lowlevel_wrapper(A->_mp_d, B->_mp_d, P->_mp_d, c_limbs, q_limbs, LOG2P/64);
	gmp_printf("C = 0x%Zx\n", A);
	mpz_set(A, C);
	gmp_montgomery_wrapper(A->_mp_d, B->_mp_d, P->_mp_d, invP->_mp_d, q_limbs, LOG2P/64);
	gmp_printf("Cphim1 = 0x%Zx\n", A);
	mpz_set(A, C);
	gmp_montgomeryCIOS_wrapper(A->_mp_d, B->_mp_d, P->_mp_d, &mip0, q_limbs, LOG2P/64);
	gmp_printf("Cphim1 = 0x%Zx\n", A);
	mpz_set(A, C);
	gmp_montgomerycustom_wrapper(A->_mp_d, B->_mp_d, P->_mp_d, invP->_mp_d, q_limbs, LOG2P/64);
	gmp_printf("Cphim1 = 0x%Zx\n\n", A);
	mpz_set(A, C);
	gmp_montgomeryparttrunc_wrapper(A->_mp_d, B->_mp_d, P->_mp_d, invP->_mp_d, q_limbs, LOG2P/64);
	gmp_printf("Cphim1 = 0x%Zx\n\n", A);*/
	
	BN_CTX *ctx = BN_CTX_new();
	BN_CTX_start(ctx);
	BIGNUM *opP = BN_CTX_get(ctx);
	BIGNUM *opA = BN_CTX_get(ctx);
	BIGNUM *opB = BN_CTX_get(ctx);
	BIGNUM *opC = BN_CTX_get(ctx);
	BN_MONT_CTX *mont_ctx = BN_MONT_CTX_new();
	
	char* stmp = mpz_get_str(NULL, 10, P);
	BN_dec2bn(&opP, stmp);
	free(stmp);
	BN_MONT_CTX_set(mont_ctx, opP, ctx);
	stmp = mpz_get_str(NULL, 10, C);
	BN_dec2bn(&opA, stmp);
	free(stmp);
	stmp = mpz_get_str(NULL, 10, B);
	BN_dec2bn(&opB, stmp);
	free(stmp);
	BN_mod_mul_montgomery(opC, opA, opB, mont_ctx, ctx);
	stmp = BN_bn2hex(opC);
	//printf("Cphim1 = 0x%s\n\n", stmp);
	free(stmp);
	
	printf("GMP montgomery: %ld\n", gmpbench(A, B, P, r, invP->_mp_d, q_limbs, W, gmp_montgomery_wrapper));
	printf("GMP montgomery CIOS:%ld\n", gmpbench(A, B, P, r, &mip0, q_limbs, W, gmp_montgomeryCIOS_wrapper));
	printf("OpenSSL montgomery CIOS:%ld\n", opensslbench(A, B, P, r, ctx, mont_ctx, W, BN_mod_mul_montgomery));
	
	mpz_clears(P, R, invP, A, B, C, NULL);
	gmp_randclear(r);
	
	BN_MONT_CTX_free(mont_ctx);
	BN_CTX_end(ctx);
	BN_CTX_free(ctx);
	
	return 0;
}

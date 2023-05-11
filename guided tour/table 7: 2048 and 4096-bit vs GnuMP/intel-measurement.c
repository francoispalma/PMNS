#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "pmns2048.h"
#include "pmns4096.h"
#include "pmns2048128.h"
#include "pmns4096128.h"
#include <gmp.h>
#include "gmp_stuff.c"

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

void do_benchgmp(uint64_t retcycles[3], const char* pstr, const uint8_t W)
{
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	mpz_t modul_p;
	mpz_init(modul_p);
	
	mpz_set_str(modul_p, pstr, 0);
	
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
	
	mpz_clears(modul_p, A, B, NULL);
	free(c_limbs);
	free(q_limbs);
	free(mip_limbs);
	gmp_randclear(r);
}

int main(void)
{
	const char *p2048 = "27281550988635662796448919634591408393040878225180401238088013120706012509944490268195902970257771551335599136987708881025959111693025021232483587405683897530724360960152562126267070784136672499918603803443651141334218499796387468846114462769965543313515054360312741752122634187095574497317619606207438286444007817917575362307092658207302609956764534108720045336517053603019146524350111615946966256259203605925378054356971762373789579733672474466221985258631945982834794491197109777205795541136241312874384112614228304693146219063565934950184249182661369862363683255030886320070261637972142236633478605252844199761123",
		*p4096 = "912269796860889774907257314466849160714212891199856033461930153027385190063358714960341003415875064237013317418693379058227236332862167370484350210518844230270038704791230532565453944927938269706678854917810025551904461314856602141494074935479078077317793050677444343496351597694127368365104399100721610624399190814780585179185630649689838734287884704088531578662737145786087309454367354511418665191985571361081027780867574148130810915279285591951586577031956385009859257036624947460607162351203755121569796831328923724787518780165875235804766106988344710692282572998046914290486942562079074088910803258544362425937497425631711239608603366098028028030386137891973660968701962407175133730274099801395096636774406971031231345716581396918871101466323932347499728769922400520092703611609413415037655339154314918456644075515882517394444818911118619726213789948570152921039618963813613108114834424072494145149368705670640697425681771434732015107724872201266346122341526620892795545198759346562843400408295301954585684932485564669674837877391709859239140735501719317619926962305353291794339493228589399485855311112216030623678529171361252699984869886501997516397292205643411658689373387601941901675598919509316760672911969406681129707254591";
	
	const uint8_t N2048 = 40, N4096 = 84;
	const uint8_t RHO2048 = 57, RHO4096 = 56;
	const uint8_t N2048128 = 18, N4096128 = 36;
	const uint8_t RHO2048128 = 118, RHO4096128 = 119;
	uint64_t cycles2048, cycles20482core, cycles4096, cycles40966core;
	uint64_t cycles2048128, cycles20481282core, cycles4096128, cycles40961284core;
	uint64_t cyclesGMP2048[3], cyclesGMP4096[3];
	
	printf("\nStarting.\nWARNING: The measures might be faulty if this is not run on a computer with 6 cores available.\n\n");
	cycles2048 = do_bench(pmns2048_montg_mult, N2048, RHO2048, 1);
	cycles2048128 = do_bench128(pmns2048128_montg_mult, N2048128, RHO2048128, 1);
	do_benchgmp(cyclesGMP2048, p2048, 1);
	cycles20482core = do_bench(pmns2048_montg_mult_2core, N2048, RHO2048, 1);
	cycles20481282core = do_bench128(pmns2048128_montg_mult_2core, N2048128, RHO2048128, 1);
	cycles4096 = do_bench(pmns4096_montg_mult, N4096, RHO4096, 1);
	cycles4096128 = do_bench128(pmns4096128_montg_mult, N4096128, RHO4096128, 1);
	do_benchgmp(cyclesGMP4096, p4096, 1);
	cycles40966core = do_bench(pmns4096_montg_mult_6core, N4096, RHO4096, 1);
	cycles40961284core = do_bench128(pmns4096128_montg_mult_4core, N4096128, RHO4096128, 1);
	
	printf("\n===============================================================================\n");
	printf("|   Method\\Size   |            2048             |            4096             |\n");
	printf("===============================================================================\n");
	printf("|        n        |             32              |             64              |\n");
	printf("===============================================================================\n");
	printf("|    Low level    |            %lu             |            %lu            |\n", cyclesGMP2048[0], cyclesGMP4096[0]);
	printf("| Classical Mont. |            %lu             |            %lu            |\n", cyclesGMP2048[1], cyclesGMP4096[1]);
	printf("|   mont. CIOS    |            %lu             |            %lu            |\n", cyclesGMP2048[2], cyclesGMP4096[2]);
	printf("===============================================================================\n");
	printf("|        n        |             %d              |             %d              |\n", N2048, N4096);
	printf("===============================================================================\n");
	printf("|      Red-64     | %lu(1-code) %s%lu(2-core)  | %lu(1-core) %lu(6-core) |\n", cycles2048, cycles20482core > 10000 ? "" : " ", cycles20482core, cycles4096, cycles40966core);
	printf("===============================================================================\n");
	printf("|     Red-128     | %lu(1-code) %lu(2-core) | %lu(1-core) %lu(4-core) |\n", cycles2048128, cycles20481282core, cycles4096128, cycles40961284core);
	printf("===============================================================================\n\n");
	return 0;
}


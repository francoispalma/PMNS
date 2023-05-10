#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "pmns1024.h"
#include "pmns2048.h"
#include "pmns4096.h"
#include "pmns8192.h"
#include "pmns8192128.h"
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
	
	const char *p1024 = "151557439827487301966960191905248420233409327748281249278640976537066266993988257219581147225805470889269899198735325674900327895538259785923482269207281672809731126255742686220820238394778878387647139546578384166978770043318856199095456814458300821186859906675966709623666377298322314372989786653793067293709",
		*p2048 = "27281550988635662796448919634591408393040878225180401238088013120706012509944490268195902970257771551335599136987708881025959111693025021232483587405683897530724360960152562126267070784136672499918603803443651141334218499796387468846114462769965543313515054360312741752122634187095574497317619606207438286444007817917575362307092658207302609956764534108720045336517053603019146524350111615946966256259203605925378054356971762373789579733672474466221985258631945982834794491197109777205795541136241312874384112614228304693146219063565934950184249182661369862363683255030886320070261637972142236633478605252844199761123",
		*p4096 ="912269796860889774907257314466849160714212891199856033461930153027385190063358714960341003415875064237013317418693379058227236332862167370484350210518844230270038704791230532565453944927938269706678854917810025551904461314856602141494074935479078077317793050677444343496351597694127368365104399100721610624399190814780585179185630649689838734287884704088531578662737145786087309454367354511418665191985571361081027780867574148130810915279285591951586577031956385009859257036624947460607162351203755121569796831328923724787518780165875235804766106988344710692282572998046914290486942562079074088910803258544362425937497425631711239608603366098028028030386137891973660968701962407175133730274099801395096636774406971031231345716581396918871101466323932347499728769922400520092703611609413415037655339154314918456644075515882517394444818911118619726213789948570152921039618963813613108114834424072494145149368705670640697425681771434732015107724872201266346122341526620892795545198759346562843400408295301954585684932485564669674837877391709859239140735501719317619926962305353291794339493228589399485855311112216030623678529171361252699984869886501997516397292205643411658689373387601941901675598919509316760672911969406681129707254591",
		*p8192 ="835787394958289162589853935177239754779053474934008317911505814396989805347791321094322353411823511640084543938464857683631256642915752103637527700570631550266519072451310075313770751512161963703580289907478240688159601071255018234281154813047216640349035233967166422505247420410975544165253962889683979900716513891874044492227291009763178239878517060415622376125594083116263932472149274146140523811576007036161139940157160561629997660341005729681013894895063780936858910932395630784075344172590667357410259149798542616081374784814477730402535540246645135395581078617752247140304244194263911515055442771281483545580422765803284720227978729380538000789258375215020176056339006311212563611054597171402619240752178172481519217435758136030378725180191933700193964271799490370348359683511576401664084395590180497269140982599939451006443449231936830650074768143734342649283929661264801499077388830403569048153155823656675542106083919975015560245918964007569496701903711706246627060487705832049436293160412402609721851533893684534479535497402461782079008631097358089031661531325713434533167892968455596745254254309282918025449619617029674750018614934058476776424561338585954217360667413917119656355631024657845277299800047809488780451814915818383271781304710484787157604565453682346210421888904126964828110369746584006150071033390197997229588403666087655239748328468656389185640456033549398537717114180006182928579192626791728869564768770877022013745455506031088967505891747566553547830273791372496416592067383206086485793981910839517766627160393689869538016550065224392636421169223571066714506801204581508640705062577807986653648287032237229192245686544058012487924279765690195980886440803821818614218364831036109702392592914929180690252914492587982543268211851699135619975798168179981367061304472058026665954953450414864073516882929830280958357753894763597994867295909833098906527801721545882400470916016340449229453758846769370051032621970146815950173434640131528668244391701151590404073995170863246727751632601280593598954529631274758189149527976922984385724994822084071141587664239995871795442948609968844207152713239603323134409143768695489087848412144162428421836583831263718113356279345989715729898695902516136117574080208289208571933285083157114678921248313946604042046826387276404653250395899720593881752463294679877259788707898450037369248936518132221616237922654160470001710037299449400903290051998137425060141724721329085977550380844165092914750509884625097469";
	
	const uint8_t N1024 = 19, N2048 = 40, N4096 = 84, N8192 = 187;
	const uint8_t RHO1024 = 58, RHO2048 = 57, RHO4096 = 56, RHO8192 = 55;
	uint64_t cycles1024, cycles2048, cycles4096, cycles8192;
	const uint8_t N8192128 = 72;
	const uint8_t RHO8192128 = 120;
	uint64_t cycles8192128;
	uint64_t cyclesGMP1024[3], cyclesGMP2048[3], cyclesGMP4096[3], cyclesGMP8192[3];
	
	printf("\nFirst checking A times B mod P for each size of PMNS\n");
	
	//1024 bit
	_mpnum __P1024__ = { .deg = 16,
		.sign = 1,
		.t = (uint64_t[]) {0xd38a4aa5f0d15c0d, 0x1a492619bacb053b, 0xa756c4c8676a2cc7, 0xf6aa607a79d82cb5, 0xcffaf035c7912dcf, 0x578088fb546e8053, 0x39f9e77fc5a7d519, 0x1d781509c1649be8, 0xe9d5984919726224, 0xdda2ff47cb3732fc, 0xfc82a0bbecb47552, 0xa681dd40133a7bee, 0xca9098d4c870e97f, 0x8e8ad494318a1464, 0x7505298b399c3183, 0xd7d330786f4f8808} },
	Gamma1024 = { .deg = 16,
		.sign = 1,
		.t = (uint64_t[]) {0x6f2f0c40b3eb32f2, 0xca5eaac5e615854b, 0x256555f023f01eb1, 0x73868edca84d52c3, 0xf56a331daeebe76b, 0x76e2ad01145cf016, 0xdae907cb49f8ad6a, 0xba4b0b174c5a96e7, 0x4ef4b080a8b91b74, 0x75759f27ca976704, 0x372684148d231ce1, 0xa5dbcc9744ca7fbc, 0x8a08d97b3e75b45e, 0xfaa82a600e4d5115, 0xcc4d0d4c581bbb13, 0x109d383e4c597848} };
	
	printf("Checking 1024 bits:\n");
	check_pmns_vs_gmp(&__P1024__, N1024, &Gamma1024, RHO1024, pmns1024_montg_mult);
	
	
	// 2048 bit
	_mpnum __P2048__ = { .deg = 32,
		.sign = 1,
		.t = (uint64_t[]) {0x6d84986d6c2a54e3, 0x7461a84e07b551e7, 0x57eede3e1105227f, 0xf30a90db56b60978, 0xdf791ec1591cf01f, 0x4d3c81d8256d7ad5, 0xfe150116172b9b76, 0x2fa279ee344d23f6, 0xd151883932c93ec3, 0x31e26b2b4e35cd7f, 0x360b972869de36f4, 0xfed99eefdbd17697, 0x146915c4e7090dec, 0x5081888c29ba9a05, 0x93dc2f0678b3d65d, 0xaa8b10a85a0fa135, 0xf1b08c7337b40680, 0x3a165ff7b7b2c9c4, 0x87765d1c3cf29e0, 0xbf34f0d73d704d44, 0x4b6483e9d7537413, 0xc9f9ee23941713f3, 0xf94b03a8ae7faa6b, 0x802a6d4d32ab0403, 0xa1df5a00f73b2141, 0x3329f2142114abfa, 0x848ddfd641a09173, 0xf1d1c68bc8533ca9, 0x458cafca01cefaab, 0x4d095596ba8388f2, 0x1817b4ad9cf09522, 0xd81c8c10e9516878} },
	Gamma2048 = { .deg = 32,
		.sign = 1,
		.t = (uint64_t[]) {0xeabc5b783fbfade1, 0xcff39cf87a3fe6ba, 0xde2d8c3d0cfe1040, 0x3d0dd05ad3550a8f, 0xc86a1cdbdd9e4b2f, 0x23c4cdff0a43060a, 0x2d2f836f0032a713, 0x67caf46d559afa9c, 0xded4ab51f7ea2fb8, 0xacf86bd7f86baa3d, 0xbf09522bcb732fe1, 0x6563e17a1985455b, 0xe7029e76626f814e, 0xcf4e0399ecd87d0b, 0x752e55aff51f45f2, 0xb1bfaeeea1d37f2f, 0x29d83e0abaa5cfd7, 0xece018f744256f8d, 0x991e51a4986de02d, 0xf19ae02986e671db, 0xfee6f47fe1fb1e20, 0x75e3de25e1a4fde, 0xeee0b961463b145a, 0x60507cee50786c68, 0xe60e8ac445bd1199, 0xa0ab51726e051c42, 0xba1da218037afc47, 0xa9c4ef8573f7b5cc, 0x193b7d1e50400c53, 0x318223c3acd4e778, 0xe74cd2401464242c, 0xa16cc487ae9b6894} };
	
	printf("Checking 2048 bits:\n");
	check_pmns_vs_gmp(&__P2048__, N2048, &Gamma2048, RHO2048, pmns2048_montg_mult);
	
	// 4096 bit
	_mpnum __P4096__ = { .deg = 64,
		.sign = 1,
		.t = (uint64_t[]) {0x146e758d65f2333f, 0xc6a7d3893b3dc80c, 0x9c52e0fdb0ec7722, 0xaddb7d55734c85ef, 0x5f2cec8c0ddb9fbe, 0x357e52a3694f6edd, 0x51131b112d51b23, 0xa1018d0a22a45dae, 0x4aaca171f8671239, 0x149580b070d596c3, 0x82ecc004eefaad55, 0x91eb02b5bd7684bb, 0x3c3b78c8ffce3ea9, 0x5140b457c79c3509, 0x8f616d7211cfc3fb, 0xfb4d3deb2c85e4c5, 0xe50ff6d59a95841e, 0xa321bf4ff1822b05, 0x28182eb1f32f559a, 0x935370913c01b32d, 0xd98becd771b420cd, 0x81e9d78e5eabf030, 0xcf674658e8046443, 0x992c35da6b27df39, 0x26ce4d56d673221c, 0x80ea6618652ff333, 0xa419909ff6a2e585, 0x3bd0879714ba98ae, 0x819a4bbe53947671, 0x764479a919023d52, 0xb182f15bf6d6559, 0x7a3189689a65e28e, 0x10d30926577bf814, 0x7881e921296b5511, 0x5a6e366541c8dd26, 0xafe5255a4ba61527, 0xd19ba91e22919473, 0x93e5ccbb6dbefd8a, 0x24f89dd9cf82b12d, 0xc52d4d5cfc668095, 0x5553bab3b85dc0ce, 0x5aba1e195636db65, 0x667eee11e5725189, 0x4259d3cc9ea8df9, 0xef97ff4ce0114116, 0x245cbb85cb3259bd, 0x4abda85b9c0993ef, 0x181bc8301d560d8e, 0xc055426486e46236, 0xa3cce3392c5a8567, 0xad3e9a28fc77edea, 0x11663397aa9bc21d, 0xbc67e671ca7f9fa7, 0xadbc5a6c39b47763, 0x369ba1e22f6dfba4, 0xe965548c7dc1e388, 0xa2845e745b51fdfd, 0x7bb9c6c4301bfcef, 0xabb70f1c294d6f3f, 0x13323906aef43c2b, 0xf93d92dd7f7f6d08, 0x9a6179fb587e7161, 0x2d04088e91323f3, 0xdf9d73ab6766c722} },
	Gamma4096 = { .deg = 64,
		.sign = 1,
		.t = (uint64_t[]) {0xb5d7fd8301e31d84, 0x878590b5e96a309f, 0xb328cffce45730d9, 0x209d45be15415af2, 0x7e30c7214c3de898, 0xef3289917e744136, 0x7b9ac1a342bee3f8, 0x963f40e3005e79b9, 0xcde04f3c26f6ffa, 0x7d110e2f7caf2b0b, 0x22a0ebd637f90020, 0x8b03873aeb3e34e7, 0x6e4936a463511524, 0x1fedb11890ac3214, 0x6e0ae4e6e8321310, 0x1f201d7a197dd220, 0x1feed647772fad84, 0x789f5171ffc17ee1, 0x47068cd865128074, 0x6fb626df85dd7032, 0x4fc66a14f0e8c31a, 0xc754703e0c88c18, 0x3a34aa67cb3c448e, 0x2c52c564acc56529, 0x61b7db18a9abdea, 0xe647de56b5f0a449, 0x785985ee8a88fd63, 0x73c1918ae6a5b782, 0x533d02f493c049e9, 0x980475087c865df5, 0xc26067e81c362ab1, 0xf5c7c3ee22910be9, 0x5be4a258a373620f, 0x6cd3477a8609f58b, 0x155e104542015ccb, 0x8d4b2623dc99e023, 0x75794ee246afc72e, 0x53d121006aed99d2, 0x96bd5d49305b5ce2, 0x9c2a463bbfb32dd2, 0x91aca706454349d7, 0x33a7390053730f92, 0x9c9a216e90e00b8a, 0x504fa70cdc7d1adf, 0x8061acf5ff0ef957, 0x1d19881a51700f47, 0x5de1e370d1ed8de6, 0xf2ba67d56f788d60, 0x12227d0d73ac4b35, 0xa9a3b0e412609f1b, 0x6334a3b63c740ca8, 0x8868856ac853678d, 0x24fcadd5574f8781, 0x9b18d075bcd47f9f, 0x9ee83f0c499420d, 0x2480eeccda36cb58, 0xef1ffed0a7b34a25, 0x53b93a445ca53818, 0x3e65c255f68a70de, 0x6d7b84fe1af77fc9, 0xfeb94df1dc3dcc68, 0x24f2005c459159e7, 0x161c9c62bc20a0e9, 0xb97c84bce9558b0c} };
	
	printf("Checking 4096 bits:\n");
	check_pmns_vs_gmp(&__P4096__, N4096, &Gamma4096, RHO4096, pmns4096_montg_mult);
	
	
	// 8192 bit
	_mpnum __P8192__ = { .deg = 128,
		.sign = 1,
		.t = (uint64_t[]) {0xc5431f03cb0792fd, 0x99d65e99ac4791c, 0x6bfa2426e7701f1a, 0x6464e653e7cbc2de, 0xc68e99232e3516d, 0x950735f1a0de93a4, 0x49eac104bb336833, 0x686d958f44df18bd, 0x75f69129f0fb179b, 0x28ee0b74aa3ab34e, 0x1d603b439ff05699, 0x8afa06a9cd6ca6a6, 0xbb72ac2fd6d73543, 0x6e4b09ade7ca6f0b, 0x3c15647063204345, 0xd0b61354ea9ef019, 0xb5b9b651b479ef09, 0xdff56ad0f82e2776, 0x9db326e8776f197b, 0x4033dbfe623bd905, 0xfb6560247c6f9685, 0x212b180de60dfe98, 0xe264db796358afb6, 0x853ce5671e375fec, 0xbb0c577d2be205b0, 0x3333a8110ac945cd, 0xc3f4fbbeab139e5b, 0x59632f79f3fc13ef, 0x3269d9464b9bd34c, 0x94fe15660c4ec41b, 0xaab3b25aaa444c88, 0x4cccdd6f1b655c9e, 0x87a512fb187dd781, 0x4756e29bbc31d57f, 0x514cda02f8904c57, 0x9a80803044a5e651, 0xcabe2c1c0757ffd0, 0xd2eb9be826d1af5e, 0xd803ddd777267984, 0xb6318bfa2dd08f89, 0xe3fdad69bde6e37f, 0x414951bf5466c2be, 0xc001f35210b7c8f4, 0x445e43a8b1955ad, 0x568f40dba6d59cdd, 0xd73ab02df64d9ff7, 0xb5fbd94c9c6ccd86, 0x1eae6c24b3c7482b, 0xe6bafb13412ce2b6, 0x596695284381a953, 0x7b1b7ec89f6c7e43, 0xf00ca5d405669517, 0x4eb4cc9e6bd03a27, 0xb177832f4ec98a7a, 0x19d3a36256b368ff, 0x92aff17deffa0137, 0x3fea6022c4325f7f, 0xf61940439b354127, 0x98cc252e76a95a5b, 0xff4c6d15c9c2640e, 0xba054a5e63ea015e, 0xa94946a928a6e207, 0xe210fe6efe7eeb3b, 0x36a3215b351ac741, 0xcfa1119a2101cf43, 0xa9182208a98ce7cc, 0x9f3d7ceedbecfdf3, 0x8ba2d7d95f97d97b, 0xd1b9a60093934305, 0xa8f61da11d814b7b, 0x1b99f9fc316c2fda, 0x442d87af138e008e, 0xfaa6b8e7e94fcf3b, 0x7100365f336db74f, 0x6fd273565cae0909, 0x342deec8c08327bc, 0xe28b64608563ece5, 0xf9664f8bd26afde6, 0x67034243c384f2ff, 0xdc1ae0f32b43954, 0xe66c3548e8df27c8, 0xfd2bb19a90fcf4b8, 0x9b2c6c9e49b6f884, 0xf4040f4f3cbe891f, 0x8488ece9123e434c, 0x959a5211ace24a76, 0xa5301ec0dc6d8436, 0xc9305c9ba2352186, 0x1b412080e6742a31, 0xd405907b95159f53, 0xe50569f69fff78fb, 0x6d3baa2f6239be29, 0xcd34685aaf7854ae, 0xa71fa93c48c37d1, 0x14d0cb11176d2f44, 0x8d1aef53b62e1a88, 0x58a09be065f522fe, 0xc082e898cd677dfb, 0x150dd50675b9e08f, 0xc93b1668079fb391, 0xa731ea5ccc46ab49, 0xfd01fb5a8ef50c97, 0x634d9b075cc143b2, 0xafbb5e56d24399d9, 0xcee0c7f5b29bfe19, 0xafbff4fc4a9b6c13, 0xb88e9a562b26ba34, 0x6f86ca54337304e1, 0xd2f975b1ee0e35ee, 0xba6166ad8bbf2af9, 0x2cccbff6893f3fb9, 0x2c17d879d51f937c, 0x5ab173aaa1b445e4, 0xde27926e82d2dfb3, 0x55a564bec12404b5, 0xfcb7a094d1a8228a, 0xd2da0bd14abd7f23, 0x4fe90779ddef8f5c, 0x26d5ab1223fd7052, 0xf842a82959a5afb7, 0x6ef789aa164abe2e, 0x8c3a7724cfa081d9, 0xbdf160b9de4ec971, 0xba834d06354a755d, 0xe2bbc2d05ad42a9c, 0x4cfd88e78569a9cd, 0xf206147857a7f11d, 0xc4290eecafd94983} },
	Gamma8192 = { .deg = 128,
		.sign = 1,
		.t = (uint64_t[]) {0xa7dc36dd3ea9db77, 0xca62d6e8b4cbd6b7, 0xae046d105ef73ec4, 0xc1313513e8cb7e2, 0x2a8e9c44b26965b4, 0x8d4cf495cd8d706f, 0xbc6561a081e5cbd9, 0x1f67222729c0d204, 0xf4f5f7a8aa52b88b, 0x820a0e8bf23268fa, 0x4d948c941f912999, 0x5b553e710def6fd0, 0x1883b49225e2c2f6, 0x5963167800615549, 0xdf90f08d901d2504, 0x2c7dc0e405d1dc, 0xa2b6650893df9191, 0x494a5a2df04362a7, 0x47f71fcd132f6e04, 0x2f0810eb6f4cbca2, 0xc92a3382fc81dc41, 0xf989020dbfb3c353, 0xff5cdde4989778c9, 0x11f9f0f30e49a1b0, 0x541e3d0f8c3fa577, 0xa7264358a4c1a723, 0xc889088e70710314, 0x875ce3bde22eb166, 0x62ada58efeef1ef, 0x132dda053b6c4ab5, 0xeb695d7c97fa88f2, 0x1c71a7ac3db06b3b, 0x54b73524a1fff68f, 0xd4d594171ef228e4, 0x2615774d5c2b0d8b, 0xaf4b34845b4fcef7, 0x3c1e96b757ec0618, 0xf13709058cf29cee, 0xe285028e169636b3, 0x5d66c0027ea77b8e, 0x758afc81bc5258bb, 0xc71c2a032bdedc6b, 0xe5898674436ce3fb, 0x875cd866de639e43, 0x8f87f1207b3ca356, 0xde44e49c2225d717, 0x50c910f96afe1aa2, 0x3fe6f138e71fbd99, 0xfecf5a7ea41957a9, 0x4cf95ed8c6c6040a, 0x80d177fa87067aec, 0x9f673f8c769ac5c4, 0xb6011d84ff3b7d00, 0x335793944bc5a4bf, 0xd50e9cb2d1c28c8c, 0x464bdb5a58de39df, 0xaecb7c92ecc0548d, 0x59dd09ac9155a5dc, 0xf21f6c80604cb2b9, 0x62a10ff73de8a6ce, 0x311d8aba20d61682, 0xab76ce8ca1e1bc0d, 0xd06d2fbf181f8d9f, 0x32be8a81de703f5c, 0x5d7feb2310b5b01c, 0x4fc3bd3db68018e8, 0x178fc86cd6e93f73, 0x3ed6475c2da178e4, 0xc91a38fec086e5e2, 0xfbe75f1396693f94, 0xe26f668b03ab2c24, 0x3ead8df206ad718f, 0x6741b3b0624caa7f, 0x8a7b852d01da93df, 0xc91e40b41e4cfebf, 0xe1fa775611430746, 0xced20da4c510740, 0x1d3a9c3920fc3e50, 0x8ab9d04a985aa4eb, 0xcb86e152626151b6, 0xfb23d896d4eefff5, 0x9344168e5123c000, 0xd6e77a98ba085eb5, 0x5aac8a27006dddfb, 0x4a4f2ede2e6a147, 0x2187097423b76326, 0x70d404d5ce45c36b, 0xe1f01df66f83f859, 0x58d37a9ee5df4239, 0xac68f594a9dfcb55, 0x443e0efc3e2228c4, 0xe36d9c29f6009be0, 0x3d1b44887a18e767, 0x6380d80841351d23, 0x90b689caba5fee1a, 0x9606db13876f80cb, 0x4f5105c4a3500d89, 0xfa2a919401abe0cc, 0x775b20b2dda6bdf0, 0x1c12ee9c815391ef, 0x63975f12838c9202, 0x143168463d77ada2, 0xf2c15156215a5f30, 0x8cd1d207ecb56288, 0x63768a4978553660, 0x22d6858c5e6d36c7, 0xa05890955c831e73, 0xe5cb0273a427cea8, 0x108e011df86758e0, 0xd9b8817a4df574da, 0xc7c09372972d20e5, 0xf4ec3bb578bd299e, 0xaaa18e38c58f0d62, 0x498fb8ed1e377432, 0x7cd453903d940014, 0x75d6dab2d9242b7b, 0xf7e28c920557d267, 0x79ffff38e9ecfc6d, 0x6d98160d6b15431, 0x49d2ce1e6d62ecaa, 0xd6a791d39891bb2a, 0x42a476d43f06ea32, 0x916eda27acd9ce74, 0xffbb567c52900730, 0xcb16acc441fad6d1, 0x36cafba5c18a26a4, 0x37a8da5bd94ea735, 0x160f8eada39fc244} };
	
	
	printf("Checking 8192 bits:\n");
	check_pmns_vs_gmp(&__P8192__, N8192, &Gamma8192, RHO8192, pmns8192_montg_mult);
	
	
	printf("\n\nNext measuring execution times for each prime size\n\n");
	
	cycles1024 = do_bench(pmns1024_montg_mult, N1024, RHO1024, 10);
	do_benchgmp(cyclesGMP1024, p1024, 10);
	printf("1024-bit done\n");
	cycles2048 = do_bench(pmns2048_montg_mult, N2048, RHO2048, 1);
	do_benchgmp(cyclesGMP2048, p2048, 1);
	printf("2048-bit done\n");
	cycles4096 = do_bench(pmns4096_montg_mult, N4096, RHO4096, 1);
	do_benchgmp(cyclesGMP4096, p4096, 1);
	printf("4096-bit done\n");
	cycles8192 = do_bench(pmns8192_montg_mult, N8192, RHO8192, 1);
	cycles8192128 = do_bench128(pmns8192128_montg_mult, N8192128, RHO8192128, 1);
	do_benchgmp(cyclesGMP8192, p8192, 1);
	
	printf("\n==================================================================\n");
	printf("|     size of p     |   1024   |   2048   |   4096   |    8192   |\n");
	printf("==================================================================\n");
	printf("|                             GnuMP                              |\n");
	printf("==================================================================\n");
	printf("|         n         |    %d    |    %d    |    %d    |    %d    |\n", 16, 32, 64, 128);
	printf("==================================================================\n");
	printf("|     Low level     |   %lu   |   %lu   |  %lu   |   %lu   |\n", cyclesGMP1024[0], cyclesGMP2048[0], cyclesGMP4096[0], cyclesGMP8192[0]);
	printf("|  Classical Mont.  |   %lu   |   %lu   |  %lu   |   %lu   |\n", cyclesGMP1024[1], cyclesGMP2048[1], cyclesGMP4096[1], cyclesGMP8192[1]);
	printf("|     mont. CIOS    |   %lu   |   %lu   |  %lu   |   %lu   |\n", cyclesGMP1024[2], cyclesGMP2048[2], cyclesGMP4096[2], cyclesGMP8192[2]);
	printf("==================================================================\n");
	printf("|                           This work                            |\n");
	printf("==================================================================\n");
	printf("|         n         |    %d    |    %d    |    %d    |    %d    |\n", N1024, N2048, N4096, N8192);
	printf("==================================================================\n");
	printf("|      Red-64       |   %lu   |   %lu   |  %lu   |  %lu   |\n", cycles1024, cycles2048, cycles4096, cycles8192);
	printf("==================================================================\n");
	printf("|         n         |     -    |    -     |    -     |    %d     |\n", N8192128);
	printf("==================================================================\n");
	printf("|      Red-128      |     -    |    -     |    -     |  %lu   |\n", cycles8192128);
	printf("==================================================================\n\n");
	return 0;
}


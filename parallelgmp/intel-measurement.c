#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "structs.h"
#include "pmns8192.h"
#include <gmp.h>
#include "gmp_stuff.c"

#define NTEST 501
#define NSAMPLES 1001

#define UNUSED(X) (void)(X)

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

void randlimbs(mp_limb_t g[])
{
	for(int i = 0; i < 128; i++)
		g[i] = randomint64();
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
			for(int soak=0; soak < W/3; soak++)
			{
				pmns_mult(c, a, b);
				pmns_mult(a, b, c);
				pmns_mult(b, c, a);
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

uint64_t do_gmpbench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[]), mp_limb_t p[], mp_limb_t mip[], const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	mp_limb_t a[128], b[128], c[128];
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randlimbs(a);
		randlimbs(b);
		gmp_mult(c, a, b, p, mip);
	}
	
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
				gmp_mult(c, a, b, p, mip);
				gmp_mult(b, c, a, p, mip);
				gmp_mult(a, b, c, p, mip);
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

mp_limb_t tmp[256];
void gmpmulmod(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	UNUSED(mip);
	mpn_mul_n(tmp, a, b, 128);
	mpn_tdiv_qr(c, a, 0, tmp, 256, p, 128);
}

mp_limb_t clolo[16][128];
mp_limb_t soak[129];
void lololohihilohihi(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	UNUSED(mip);
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i < 4; i++)
		mpn_mul_n(clolo[i], a + 64 *((i & 2) > 0), b + 64 * (i & 1), 64);
	tmp[192] += mpn_add_n(tmp + 64, clolo[1], clolo[2], 128);
	tmp[128] += mpn_add_n(tmp, tmp, clolo[0], 128);
	/*tmp[192] += mpn_add_n(tmp + 64, tmp + 64, clolo[1], 128);
	tmp[192] += mpn_add_n(tmp + 64, tmp + 64, clolo[2], 128);*/
	mpn_add_n(tmp + 128, tmp + 128, clolo[3], 128);
	mpn_tdiv_qr(soak, c, 0, tmp, 256, p, 128);
}

const uint8_t ctabcorr[] = {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6},
	atabcorr[] = {0, 0, 1, 2, 1, 0, 0, 1, 2, 3, 3, 2, 1, 2, 3, 3},
	btabcorr[] = {0, 1, 0, 0, 1, 2, 3, 2, 1, 0, 1, 2, 3, 3, 2, 3};

void eightsplit(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	UNUSED(mip);
	#pragma omp parallel for num_threads(8)
	for(int i = 0; i < 8; i++)
	{
		mpn_mul_n(clolo[2*i], a + 32 * atabcorr[2*i], b + 32 * btabcorr[2*i], 32);
		mpn_mul_n(clolo[2*i+1], a + 32 * atabcorr[2*i+1], b + 32 * btabcorr[2*i+1], 32);
	}
	
	mpn_add_n(tmp, tmp, clolo[0], 64);
	for(int i = 1; i < 15; i++)
	{
		tmp[64 + ctabcorr[i]*32] += mpn_add_n(tmp + ctabcorr[i]*32, tmp + ctabcorr[i]*32, clolo[i], 64);
	}
	mpn_add_n(tmp + 192, tmp + 192, clolo[15], 64);
	mpn_tdiv_qr(soak, c, 0, tmp, 256, p, 128);
}

mp_limb_t cpad[4][256];
void cutred(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	UNUSED(mip);
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i < 4; i++)
	{
		mpn_mul_n(cpad[i] + (((i & 2) > 0) + (i&1)) * 64, a + 64 *((i & 2) > 0), b + 64 * (i & 1), 64);
		mpn_tdiv_qr(soak, clolo[i], 0, cpad[i], 256, p, 128);
	}
	tmp[129] = mpn_add_n(c, clolo[0], clolo[1], 128);
	tmp[129] += mpn_add_n(c, c, clolo[2], 128);
	tmp[129] += mpn_add_n(c, c, clolo[3], 128);
	while(tmp[129] > 0)
	{
		tmp[129] -= mpn_sub_n(c, c, p, 128);
	}
	
	/*tmp[192] += mpn_add_n(tmp + 64, clolo[1], clolo[2], 128);
	tmp[128] += mpn_add_n(tmp, tmp, clolo[0], 128);
	mpn_add_n(tmp + 128, tmp + 128, clolo[3], 128);
	mpn_tdiv_qr(soak, c, 0, tmp, 256, p, 128);*/
}

mp_limb_t aux1[3][128], aux2[128], tmp2[256];
void splitmontg(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i < 4; i++)
		mpn_mul_n(clolo[i], a + 64 *((i & 2) > 0), b + 64 * (i & 1), 64);
	tmp[192] += mpn_add_n(tmp + 64, clolo[1], clolo[2], 128);
	tmp[128] += mpn_add_n(tmp, tmp, clolo[0], 128);
	mpn_add_n(tmp + 128, tmp + 128, clolo[3], 128);
	
	// tmp = a * b
	
	mpz_t display;
	display->_mp_alloc = 256;
	//display->_mp_size = 256;
	//display->_mp_d = tmp;
	
	//gmp_printf("\n0x%Zx\n", display);
	
	//mpn_mul_n(tmp, a, b, 128);
	
	//display->_mp_d = tmp;
	
	//gmp_printf("\n0x%Zx\n", display);
	
	//exit(1);
	
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i < 3; i++)
	{
		if(i == 0)
			mpn_mul_n(aux1[0], tmp, mip, 64);
		if(i == 1)
			mpn_mullo_n(aux1[1], tmp + 64, mip, 64);
		if(i == 2)
			mpn_mullo_n(aux1[2], tmp, mip + 64, 64);
	}
	mpn_add_n(aux2 + 64, aux1[1], aux1[2], 64);
	mpn_add_n(aux2, aux2, aux1[0], 128);
	
	/*display->_mp_size = 128;
	display->_mp_d = aux2;
	
	gmp_printf("\n0x%Zx\n", display);
	
	mpn_mullo_n(aux2, tmp, mip, 128);
	
	display->_mp_d = aux2;
	
	gmp_printf("\n0x%Zx\n", display);
	
	exit(1);*/
	
	// aux2 = tmp * p^-1 mod 2^128
	
	/*#pragma omp parallel for num_threads(4)
	for(int i = 1; i < 4; i++)
		mpn_mul_n(clolo[i], aux2 + 64 *((i & 2) > 0), p + 64 * (i & 1), 64);
	c[64] += mpn_add_n(c, clolo[1] + 64, clolo[2] + 64, 64) + 1;
	mpn_add_n(c, c, clolo[3], 128);
	mpn_add_n(c, c, tmp + 128, 128);
	
	display->_mp_size = 128;
	display->_mp_d = c;
	
	gmp_printf("\n0x%Zx\n", display);
	
	mpn_mul_n(tmp, aux2, p, 128);
	
	display->_mp_d = tmp + 128;
	
	gmp_printf("\n0x%Zx\n", display);
	
	exit(1);*/
	//mpn_mul_n(tmp, aux2, p, 128);
	
	#pragma omp parallel for num_threads(4)
	for(int i = 1; i < 4; i++)
		mpn_mul_n(clolo[i], aux2 + 64 *((i & 2) > 0), p + 64 * (i & 1), 64);
	tmp2[192] += mpn_add_n(tmp2 + 64, clolo[1], clolo[2], 128);
	//tmp2[128] += mpn_add_n(tmp2, tmp2, clolo[0], 128);
	clolo[3][0] += 1;
	mpn_add_n(tmp2 + 128, tmp2 + 128, clolo[3], 128);
	
	mpn_add_n(c, tmp + 128, tmp2 + 128, 128);
}

void gmpmontgom(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	mpn_mont_mul_red_n(c, a, b, p, mip, 128);
}

void gmpmontgomCIOS(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], mp_limb_t p[], mp_limb_t mip[])
{
	mpn_mont_mul_red_1(c, a, b, p, *mip, 128);
}

void do_benchgmp(uint64_t retcycles[5], const char* pstr, const uint8_t W)
{
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	mpz_t modul_p;
	mpz_init(modul_p);
	
	mpz_set_str(modul_p, pstr, 0);
	
	int nb_limbs = mpz_size(modul_p);
	
	mp_limb_t *p_limbs, *c_limbs, *mip_limbs;
	
	p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	mpn_binvert(mip_limbs, p_limbs, nb_limbs, c_limbs);
	
	mp_limb_t a[128], b[128], c[128];
	
	randlimbs(a);
	randlimbs(b);
	//(c, a, b, p_limbs, mip_limbs);
	
	mpz_t display;
	display->_mp_alloc = 128;
	display->_mp_size = 128;
	display->_mp_d = c;
	
	//gmp_printf("\n0x%Zx\n", display);
	
	//splitmontg(c, a, b, p_limbs, mip_limbs);
	
	display->_mp_d = c;
	
	//gmp_printf("0x%Zx\n\n", display);
	
	retcycles[0] = do_gmpbench(gmpmulmod, p_limbs, NULL, W);
	
	printf("low level gmp done\n");
	
	retcycles[3] = do_gmpbench(lololohihilohihi, p_limbs, NULL, W);
	
	printf("lololohihilohihi done\n");
	
	/*retcycles[4] = do_gmpbench(cutred, p_limbs, NULL, W);
	
	printf("cutred done\n");*/
	
	/*retcycles[4] = do_gmpbench(eightsplit, p_limbs, NULL, W);
	
	printf("16-split done\n");*/
	
	
	
	retcycles[4] = do_gmpbench(splitmontg, p_limbs, mip_limbs, W);
	
	printf("splitmontg done\n");
	
	retcycles[1] = do_gmpbench(gmpmontgom, p_limbs, mip_limbs, W);
	
	printf("montgomery gmp done\n");
	
	mp_limb_t mip0;
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	retcycles[2] = do_gmpbench(gmpmontgomCIOS, p_limbs, &mip0, W);
	
	printf("montgomeryCIOS gmp done\n");
	
	mpz_clear(modul_p);
	free(mip_limbs);
	free(c_limbs);
	gmp_randclear(r);
}

int main(void)
{
	const char *p8192 ="835787394958289162589853935177239754779053474934008317911505814396989805347791321094322353411823511640084543938464857683631256642915752103637527700570631550266519072451310075313770751512161963703580289907478240688159601071255018234281154813047216640349035233967166422505247420410975544165253962889683979900716513891874044492227291009763178239878517060415622376125594083116263932472149274146140523811576007036161139940157160561629997660341005729681013894895063780936858910932395630784075344172590667357410259149798542616081374784814477730402535540246645135395581078617752247140304244194263911515055442771281483545580422765803284720227978729380538000789258375215020176056339006311212563611054597171402619240752178172481519217435758136030378725180191933700193964271799490370348359683511576401664084395590180497269140982599939451006443449231936830650074768143734342649283929661264801499077388830403569048153155823656675542106083919975015560245918964007569496701903711706246627060487705832049436293160412402609721851533893684534479535497402461782079008631097358089031661531325713434533167892968455596745254254309282918025449619617029674750018614934058476776424561338585954217360667413917119656355631024657845277299800047809488780451814915818383271781304710484787157604565453682346210421888904126964828110369746584006150071033390197997229588403666087655239748328468656389185640456033549398537717114180006182928579192626791728869564768770877022013745455506031088967505891747566553547830273791372496416592067383206086485793981910839517766627160393689869538016550065224392636421169223571066714506801204581508640705062577807986653648287032237229192245686544058012487924279765690195980886440803821818614218364831036109702392592914929180690252914492587982543268211851699135619975798168179981367061304472058026665954953450414864073516882929830280958357753894763597994867295909833098906527801721545882400470916016340449229453758846769370051032621970146815950173434640131528668244391701151590404073995170863246727751632601280593598954529631274758189149527976922984385724994822084071141587664239995871795442948609968844207152713239603323134409143768695489087848412144162428421836583831263718113356279345989715729898695902516136117574080208289208571933285083157114678921248313946604042046826387276404653250395899720593881752463294679877259788707898450037369248936518132221616237922654160470001710037299449400903290051998137425060141724721329085977550380844165092914750509884625097469";
	
	const uint8_t N8192 = 187;
	const uint8_t RHO8192 = 55;
	uint64_t cycles81925core, cycles81926core, cycles81928core;
	uint64_t cyclesGMP8192[5];
	
	UNUSED(cycles81925core); UNUSED(cycles81926core);
	
	printf("\nStarting.\nWARNING: The measures might be faulty if this is not run on a computer with 8 cores available.\n\n");
	/*cycles81925core = do_bench(pmns8192_montg_mult_5core, N8192, RHO8192, 3);
	printf("5-core done\n");
	cycles81926core = do_bench(pmns8192_montg_mult_6core, N8192, RHO8192, 3);
	printf("6-core done\n");*/
	cycles81928core = do_bench(pmns8192_montg_mult_8core, N8192, RHO8192, 3);
	printf("PMNS done\n");
	do_benchgmp(cyclesGMP8192, p8192, 3);
	
	printf("\n============================================================================\n");
	printf("|  PMNS  | Low level | Classical Mont. | Mont. CIOS | 4 core low | cut red |\n");
	printf("=============================================================================\n");
	printf("| %lu  |   %lu   |      %lu      |   %lu    |   %lu   |   %lu   |\n", cycles81928core, cyclesGMP8192[0], cyclesGMP8192[1], cyclesGMP8192[2], cyclesGMP8192[3], cyclesGMP8192[4]);
	printf("=============================================================================\n\n");
	return 0;
}


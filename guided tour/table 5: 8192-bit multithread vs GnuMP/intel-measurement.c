#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

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
	const char *p8192 ="835787394958289162589853935177239754779053474934008317911505814396989805347791321094322353411823511640084543938464857683631256642915752103637527700570631550266519072451310075313770751512161963703580289907478240688159601071255018234281154813047216640349035233967166422505247420410975544165253962889683979900716513891874044492227291009763178239878517060415622376125594083116263932472149274146140523811576007036161139940157160561629997660341005729681013894895063780936858910932395630784075344172590667357410259149798542616081374784814477730402535540246645135395581078617752247140304244194263911515055442771281483545580422765803284720227978729380538000789258375215020176056339006311212563611054597171402619240752178172481519217435758136030378725180191933700193964271799490370348359683511576401664084395590180497269140982599939451006443449231936830650074768143734342649283929661264801499077388830403569048153155823656675542106083919975015560245918964007569496701903711706246627060487705832049436293160412402609721851533893684534479535497402461782079008631097358089031661531325713434533167892968455596745254254309282918025449619617029674750018614934058476776424561338585954217360667413917119656355631024657845277299800047809488780451814915818383271781304710484787157604565453682346210421888904126964828110369746584006150071033390197997229588403666087655239748328468656389185640456033549398537717114180006182928579192626791728869564768770877022013745455506031088967505891747566553547830273791372496416592067383206086485793981910839517766627160393689869538016550065224392636421169223571066714506801204581508640705062577807986653648287032237229192245686544058012487924279765690195980886440803821818614218364831036109702392592914929180690252914492587982543268211851699135619975798168179981367061304472058026665954953450414864073516882929830280958357753894763597994867295909833098906527801721545882400470916016340449229453758846769370051032621970146815950173434640131528668244391701151590404073995170863246727751632601280593598954529631274758189149527976922984385724994822084071141587664239995871795442948609968844207152713239603323134409143768695489087848412144162428421836583831263718113356279345989715729898695902516136117574080208289208571933285083157114678921248313946604042046826387276404653250395899720593881752463294679877259788707898450037369248936518132221616237922654160470001710037299449400903290051998137425060141724721329085977550380844165092914750509884625097469";
	
	const uint8_t N8192 = 187;
	const uint8_t RHO8192 = 55;
	uint64_t cycles81925core, cycles81926core, cycles81928core;
	uint64_t cyclesGMP8192[3];
	
	printf("\nStarting.\nWARNING: The measures might be faulty if this is not run on a computer with 8 cores available.\n\n");
	cycles81925core = do_bench(pmns8192_montg_mult_5core, N8192, RHO8192, 1);
	printf("5-core done\n");
	cycles81926core = do_bench(pmns8192_montg_mult_6core, N8192, RHO8192, 1);
	printf("6-core done\n");
	cycles81928core = do_bench(pmns8192_montg_mult_8core, N8192, RHO8192, 1);
	printf("8-core done\n");
	do_benchgmp(cyclesGMP8192, p8192, 1);
	
	printf("\n=======================================================================\n");
	printf("| 5-core | 6-core | 8-core | Low level | Classical Mont. | Mont. CIOS |\n");
	printf("=======================================================================\n");
	printf("|          n=187           |                  n=128                   |\n");
	printf("=======================================================================\n");
	printf("| %lu  | %lu  | %lu  |   %lu   |      %lu      |   %lu    |\n", cycles81925core, cycles81926core, cycles81928core, cyclesGMP8192[0], cyclesGMP8192[1], cyclesGMP8192[2]);
	printf("=======================================================================\n\n");
	return 0;
}


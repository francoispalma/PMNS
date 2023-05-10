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
	const char *p8192 ="1070058483015768003992514995681407384736105245325727180401608020840546719045322896670274718579753308994132592352036564047889209718860972559848058359905359070743952082947838302574022876520450636974049003463301704195697960709765141039568495293695136192758999339502127133882319941529451796297044849595644437971236863428443363680805491362266433712675621460184320617659320460107526260463969037266293074126575532864102981488441490057596898820706304612102110837173519985040496304742454378491616432567058769140725710073774879320474157325243985505544325028549782831857249784087539468167602044329252077348196978247148961428721371421963159879460499891334368177797537787541926116036688836695948338800106337532819156934624259783840105344092308693103117621016140094268703125458353260482416828545951808513108547849191137027371582442720268511113826041987643074906423574300067152315006841765846308554202747828141295254849237556387377711885149696584856600988532812292392234357808367641086137625025109535003287220090833745013438023394862922204887632492952107047717070107708049704271537294620833734515961184826445806847767417220752598417006839515864079429961063585205172485445168888419544163172980625495403396185560717355518647797710123041617889392779999986722769008233319791396809091379853401526927694854084619306368078899669843190620102807806106038825335722315297655938924410065008952487367557289542183461874132260314883579690750018678983233417598081834024230491828338996100182120355123125747113697904626233224179379850760936774166642976205457602275999773563757039255356468286891278903201727397937169829260556937590010249356342626850734218479955354576352143095203081363122396303303080663507639846827574515764834427907331626675550269418914420741103732346649209554017731666388977925921015649887875434343292469078816802186047538908754103232122010626741650197294530588025672576777904486800872611601654514158730329294768302786645269746870776841865698137878165773226110730608163616601299030102047995490720208315541881343385926533422024953323940227808465617350965530688342249105370261297027270326443930449507281007689502916874977390190133154703287047109746245707959234603073904229979479212421873111576990291362203282682224161149815073611444819337858079884767524630453949485486111155406554722075732638710215373003073367808581430466272862300818033911371863392358332046911489793971578804180505962399163386370720641691529319513004479900146636709616644715426057695152225842830786529382504260109713";
	
	const uint8_t N8192 = 189;
	const uint8_t RHO8192 = 55;
	uint64_t cycles8192;
	uint64_t cyclesGMP8192[3];
	
	_mpnum __P8192__ = { .deg = 128,
		.sign = 1,
		.t = (uint64_t[]) {0x5807e79c62830591, 0x61c742922aa3de7c, 0xbe99e94ba74e0f67, 0xb7f4db341dfa907a, 0x17d031844f417531, 0xc4b9925f6e8d55d1, 0x33c2efd940261b6b, 0x37494d38b0ed13f, 0xbca1493873390247, 0xd39c62d16ca57b74, 0xf1346261804fab64, 0x5103df3612b0b2a1, 0x3e69ee1c4183f202, 0x58b1b97c39e53f5c, 0x3582864552b75a16, 0x6eb4ac30d5db2e8e, 0x2080b5e135a5b146, 0x673e8cc25093a581, 0x9ef8ea420677af5a, 0x928242c40caa94da, 0xbfb64c5574426b22, 0xeab5de759bd55919, 0x4ffe4e3f26b424f9, 0x359b7a8c0c94a5fc, 0x98364926491c343, 0xf2a32342bffeecd3, 0xc0818cb4b7fa7c50, 0x9104814ef1d8ce4d, 0x3594a5a25c7ebcfd, 0x26e62da391286810, 0x39277ca0ac69dc35, 0xc6550cc10b244328, 0x834f818934c36ee2, 0xb0f853fb11000c94, 0xa4fc1827dbb1dd18, 0xcf51a6d2d52809f7, 0x7c60c7d8f3d17db1, 0xcad533f2c3e70388, 0x858d29da11f04009, 0xea995034a8539346, 0x3f226bd8709a1366, 0x2f9fe5a4b1af65ef, 0x33b658d406c7c82e, 0xf30fa36b23126b7, 0xf780790576ccf0a1, 0x5b268f3f6f2477a3, 0x298963eb614e4d2e, 0x490c76b90591a42d, 0xaff28f82eef9232f, 0x55604a98a1c8f582, 0x2109d7470a3467c6, 0x2298f7ed111c6757, 0x849c4fa695eef154, 0x3d36d977c189292b, 0xfb04a5929a59e9e, 0xea764b86f63e783a, 0xa817dca21ff7b284, 0x3c273449b672d4e6, 0xb34397bf1e57d75f, 0x3bc6f71476f44fce, 0xd2a132c951f102bd, 0xf095745a521211fc, 0x6c0024e1408d3d7a, 0xd652c26776024569, 0xcbffd3c55ccc4e26, 0x2c48553547cd4316, 0xc46d688f65d45691, 0x8b1bbcbfce2ec4bc, 0xf6fc94ceef0aa908, 0x9a4f5d3ab6df520f, 0xea6f6464f8737dce, 0x9991a930aded985a, 0xa7081b76afa3daec, 0x8cd3ca330b5f5da1, 0x2129c1ee04205dec, 0x82b63dc77969e969, 0x253c77200978ce80, 0x9f2cb13f7cd4e689, 0x5d6115adf56150e2, 0xc8f1e0670473bbff, 0xc9540f005482efdb, 0x399400dd5157f954, 0xf22f6bfb4807b9c2, 0xc43d8c480a817ef9, 0x48b2f42c214bbefa, 0xbc70102be09889b9, 0x69024774a6986de6, 0x66e9ef2bc4bb3683, 0xec6caf7da460e6bb, 0x50d106c51b419db3, 0xb25be2a3a388e32b, 0xa5d1340be41a933b, 0xecaecf7784a50780, 0xd725db32e6d26b62, 0x1535ac1ce06d076a, 0xd91633d589ff4995, 0x335d415a1d307a9a, 0x58c3ae75d5ba0dc2, 0xfe448b2bd7ef4db9, 0xc9c8f079d59d086a, 0xe11288049858134d, 0xf3dc4f0865189053, 0xa8f4f5050d847ba1, 0x2971ac5396a27c63, 0xfce73d147f0756b9, 0x9ba0c214aefb738c, 0xaaaa93e03612101c, 0x63779aad51191a1, 0x1213eb14ba1d089f, 0x9ffe8cf39e2163af, 0x1e7d12216138eb4, 0x1e1e4b0c81c00f91, 0xa462d27e78b462ae, 0xd5e3d50063defba3, 0x99f25a24e183eba4, 0xc72149270e663485, 0x3aa6136219d7292, 0x422e6666e420fa10, 0x28c729f2c2521b60, 0xdd33e3c40561a29c, 0x70811b83f7ba10c5, 0x499ed6483c9100b8, 0x8c57d8b9b23d4236, 0x17b1b5e141e292ef, 0x88e336286fa971ad, 0x3561b8fccc1bea42, 0x5fa4429dd4c9207d, 0xfb24e481c9a110c7} },
	Gamma8192 = { .deg = 128,
		.sign = 1,
		.t = (uint64_t[]) {0xb84cd9544966a5e8, 0xda68e9184bba8426, 0xefed98760533450, 0xab226a2cfc5aa90, 0x6d956e28926f49a4, 0x93b2f6629b215fdb, 0x4e029b0a183d377, 0xae7a0368fa5e7db4, 0x9e8f107f04d39629, 0xcbcff7351f516cfa, 0x2a1fc8b5bf7a8317, 0xd520ccb845b2b573, 0x2994735279782092, 0x774fdae5e619cbf5, 0xda358bdc9d6fb3e7, 0xf40f03c48389b5a2, 0x663b65c22f1dd802, 0xeb18360d347f11a9, 0x7b0916cb3af9d956, 0x1f52cab4102fe259, 0xe5f674b11e6cfc96, 0xda87f142d4d9c0de, 0xc0c727aa18390e49, 0x9a6d8bb66a2f1c5b, 0xd3b559edc267222a, 0x2bc906e0a1d05b28, 0x742fe2ba9c91f1f1, 0x94b7ee182492e2d1, 0x18ceeab610863502, 0x48fd568e2e1b8d0c, 0xbd68ab1488e46f7e, 0x274495d33ad7b850, 0x6c377b1f1dd8355b, 0xae1ea8442b2bac5a, 0x9cb70b52c74a2aa2, 0xcadb8864bc32beed, 0xe766e6384ff35e74, 0xeac6439474128655, 0x743b16fc1b2b8f7d, 0x16145e23c65f11d0, 0x3304b2dc46d6dc9c, 0x23edc16ea29cf26e, 0x76c6d8c96fbb58be, 0x4e55c20bc0051cf8, 0x2678f1a37f451e3, 0x1b3f09e1405f8c6e, 0xc4e0b162fe6cab0b, 0x36cc8555e06dde44, 0x7ab47bea8e93f0bd, 0x1c9f2e24bfbe1dc7, 0xe18f8b674b581522, 0x5593ba3ee97a790c, 0xd51f186582cfe304, 0xddad4951e7fdb3b0, 0x744c8cf33901384f, 0x631d1d837777df88, 0xfec5a4bad12bbf3f, 0x5343b92c8175586a, 0x3a9283f703b3e6a0, 0xc246aa79b1b3fce7, 0xba98cb2ac2e21ccd, 0x4434dfd905ff461a, 0x8a099ab191cf8e7f, 0x45bd57614e942431, 0x89adb974afa37496, 0x4984a649b6679d4d, 0x751985b56d158a01, 0x19b08c499b438f22, 0x718a3b153c3b8cbb, 0xcf82097618392bd2, 0x7ed868acf1af6be2, 0xa081191e07df45bb, 0x938c6f5e3a7980f9, 0x296493a297a1eaa4, 0xfa2830ce7240ea6f, 0x301440058d7cd03d, 0x7676f29a63b37bc6, 0xa7836d03fc3aa16c, 0xa8eb9ef33a40eb82, 0xc7ae8b422eca63d4, 0x17a1815485c84071, 0x4e054347c5496b07, 0xe0fae410abe3debf, 0xf9d85322f0be44fc, 0x7e1deafcc994eaf4, 0xfcbbe979e681d5ce, 0xce60be1096579526, 0x8654f9e37e2a8d9a, 0x64277fdfd02da9f2, 0x14cc4a017b3c1724, 0x3ca481e741d3c2e3, 0x6d2cde1ef28d01f1, 0xcdda6d1722aebb0d, 0xfd2767ef3899b239, 0x294b50c645f72208, 0xd72819703228a1ed, 0x3282fa2ccf1a6a3a, 0x173974ee291be584, 0x128b209d1692fbd2, 0xc0e950af570cef71, 0x9aa8913d9ec3973e, 0x3f55e20410f17ccc, 0x8f2182d9588c522a, 0x6cd4ed0103284681, 0xdc3bd2ef3d508af, 0x192f376e456dd378, 0xc3bb0bf4b51ac882, 0x994b55db6381e46e, 0x44132d0e31856573, 0x18e958f24d89e788, 0x704cb7dbeb6376da, 0xf4f3e2da2448867f, 0xe00e359aa9141dda, 0xe947e900422b65a9, 0x6c29d43a8601d982, 0xc3da0aac392519cc, 0x402c87c79925124c, 0x90b6f3bddefaeb2d, 0x6b1f46aca86ed2e4, 0x37e80796863acd6c, 0x685af5961c9cff37, 0xdcb6e7788e67c1da, 0x1c627849050d315, 0xf650a650f8796b7f, 0x206a5ddfa0ee7720, 0xfacb555504ac963b, 0x178ae4446950c27f, 0x7caf77402db035b8} };
	
	
	printf("\nFirst showing that the results are correct\n");
	check_pmns_vs_gmp(&__P8192__, N8192, &Gamma8192, RHO8192, pmns8192_montg_mult);
	
	printf("\nStarting.\nWARNING: The measures might be faulty if this is not run on a computer with 6 cores available.\n\n");
	cycles8192 = do_bench(pmns8192_montg_mult, N8192, RHO8192, 1);
	do_benchgmp(cyclesGMP8192, p8192, 1);
	
	printf("\n\n=====================================================\n");
	printf("| 6-core | Low level | Classical Mont. | Mont. CIOS |\n");
	printf("=====================================================\n");
	printf("| n=189  |                  n=128                   |\n");
	printf("=====================================================\n");
	printf("| %lu  |   %lu   |      %lu      |   %lu    |\n", cycles8192, cyclesGMP8192[0], cyclesGMP8192[1], cyclesGMP8192[2]);
	printf("=====================================================\n\n");
	return 0;
}


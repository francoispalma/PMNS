#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <immintrin.h>
#include "pmns.h"
#include "params.h"

#define NTEST 501
#define NSAMPLES 1001

/*typedef union fvector {
	float f[8];
	__m256 m;
} fvector;*/

#if N == 5

//256bits
const int_fast64_t translaThi[5] = {-34513430600178, -158296399171563, 271236604895078, 469474936668990, 235699262665207};
const uint_fast64_t translaTlo[5] = {11119762830108874752u, 9566583731167263744u, 6433920378123599872u, 14573075807047974912u, 11679486863747164160u};
const float floL1[5][5] = { {-4.01698188716961e-16, 1.11580929471261e-16, -1.41508995162515e-16, 2.07247224924139e-16, -5.10240321873274e-17}, {-6.38176892177692e-18, 2.60189193554447e-16, 9.56662954528780e-17, -9.04849629751872e-17, 2.07247224924139e-16}, { 6.38202933186450e-17, 1.02048064374655e-16, -3.50674156529634e-16, -1.11580929471261e-16, -9.04849629751872e-17}, {-3.94598919796291e-16, -4.14494449848279e-16, -2.13628993845916e-16, -2.60189193554447e-16, -1.11580929471261e-16}, {-6.88566026486906e-17, 1.80969925950374e-16, 1.54305256293832e-16, -1.02048064374655e-16, -2.60189193554447e-16} };
const int_fast32_t ratL1[5][5] = {{-7587865, 2107705, -2673030, 3914790, -963817}, {-120549, 4914835, 1807086, -1709213, 3914790}, {1205532, 1927634, -6624048, -2107705, -1709213}, {-7453763, -7829579, -4035338, -4914835, -2107705}, {-1300665, 3418426, 2914744, -1927634, -4914835}};
/*const fvector fvectL1[5] = {{-4.01698188716961e-16, -6.38176892177692e-18, 6.38202933186450e-17, -3.94598919796291e-16, -6.88566026486906e-17, 0, 0, 0}, {1.11580929471261e-16, 2.60189193554447e-16, 1.02048064374655e-16, -4.14494449848279e-16, 1.80969925950374e-16, 0, 0, 0}, {-1.41508995162515e-16, 9.56662954528780e-17, -3.50674156529634e-16, -2.13628993845916e-16, 1.54305256293832e-16, 0, 0, 0}, {2.07247224924139e-16, -9.04849629751872e-17, -1.11580929471261e-16, -2.60189193554447e-16, -1.02048064374655e-16, 0, 0, 0}, {-5.10240321873274e-17, 2.07247224924139e-16, -9.04849629751872e-17, -1.11580929471261e-16, -2.60189193554447e-16, 0, 0, 0}};*/


#elif N == 10

//512bits
const int_fast64_t translaThi[10] = { 1950646434596335, -2358957112726445, 781801709903168, -451066099939832, 339331992545606, -1717186922101434, 1741450998587110, 374609625217838, 745474093595848, -86221132017849 };
const uint_fast64_t translaTlo[10] = {2012500658081156096u, 8421606028041383936u, 12527811339647841280u, 9336941486971077632u, 8373280123021256704u, 13289935759032327168u, 7796273899372325888u, 16693204049766080512u, 15910050057615589376u, 14041682279933967360u};
const float floL1[10][10] = {{-2.60235926079914e-16, -4.48135474347700e-17, -1.86592740811087e-16, -2.42826940139049e-16, 1.83172939478171e-16, -4.85568192946713e-17, 1.65559096009531e-16, -2.08372904023109e-16, -1.65610447379017e-16, 1.61808868073163e-16}, {-6.61193113938249e-17, 2.60235926079914e-16, -2.33568547497327e-17, -6.61706627633112e-17, 9.33703667294414e-17, -2.17801632120228e-17, 3.73983351402127e-16, -6.99208906996792e-17, 3.44540361159399e-17, 4.85568192946713e-17}, {-1.49604566057404e-16, 6.61193113938249e-17, -4.52296392417850e-17, 2.58832821460662e-16, -2.38455762867891e-16, -4.65640359499465e-17, 3.54668545837393e-17, -6.65937106467928e-17, -3.75022793636804e-18, 2.17801632120228e-17}, {1.95641954854312e-16, 1.49604566057404e-16, 2.58485437564737e-16, 2.27358581501684e-16, -1.95552754438784e-17, -2.13640714050079e-17, 7.03439385831609e-17, 2.13671890129967e-16, -3.25426532107455e-16, 4.65640359499465e-17}, {3.23617736146326e-16, -1.95641954854312e-16, -2.15480686091097e-16, 6.85351426220314e-17, -1.28240494652396e-16, -4.48135474347700e-17, 1.11754641977488e-16, 4.47552399888170e-17, -1.36866913717165e-17, 2.13640714050079e-17}, {9.71136385893427e-17, -3.23617736146326e-16, 3.86717072288091e-17, 1.95181589195114e-16, 2.40455502289082e-16, 2.60235926079914e-16, -3.10685486171005e-17, 1.04791018622634e-16, -2.37799026332144e-17, 4.48135474347700e-17}, {4.35603264240456e-17, -9.71136385893427e-17, -8.50105948318028e-17, -1.12881248262693e-17, 6.33818100664122e-17, 6.61193113938249e-17, -8.10111159894196e-17, 6.45939712256012e-17, -9.03905705724802e-17, -2.60235926079914e-16}, {9.31280718998929e-17, -4.35603264240456e-17, -6.18564698981885e-17, -7.82736146620069e-17, 3.09943271955178e-17, 1.49604566057404e-16, 2.57965993468790e-17, -2.57498424752501e-16, 7.58820960518705e-17, -6.61193113938249e-17}, {4.27281428100157e-17, -9.31280718998929e-17, 3.76108663614387e-16, 1.44406838208765e-16, -1.06044239633358e-16, -1.95641954854312e-16, 1.81616328700630e-16, 5.24909274680614e-17, -1.79224810090494e-16, -1.49604566057404e-16}, {8.96270948695401e-17, -4.27281428100157e-17, -1.42088642689015e-16, 9.20186134796765e-17, 2.88770026754205e-16, -3.23617736146326e-16, 1.26733882622433e-16, -2.39202281278358e-16, -9.19159107407039e-17, 1.95641954854312e-16}};
const int_fast32_t ratL1[10][10] = {{-4915718, -846504, -3524638, -4586872, 3460039, -917213, 3127323, -3936053, -3128293, 3056484}, {-1248959, 4915718, -441199, -1249929, 1763717, -411416, 7064346, -1320769, 650819, 917213}, {-2825951, 1248959, -854364, 4889214, -4504303, -879570, 669950, -1257920, -70840, 411416}, {3695573, 2825951, 4882652, 4294683, -369389, -403556, 1328760, 4036148, -6147134, 879570}, {6112967, -3695573, -4070316, 1294593, -2422395, -846504, 2110986, 845403, -258535, 403556}, {1834425, -6112967, 730488, 3686876, 4542077, 4915718, -586869, 1979447, -449190, 846504}, {822832, -1834425, -1605805, -213227, 1197249, 1248959, -1530257, 1220146, -1707430, -4915718}, {1759140, -822832, -1168436, -1478547, 585467, 2825951, 487284, -4864008, 1433373, -1248959}, {807112, -1759140, 7104492, 2727769, -2003120, -3695573, 3430636, 991526, -3385461, -2825951}, {1693008, -807112, -2683979, 1738183, 5454712, -6112967, 2393936, -4518404, -1736243, 3695573}};

#endif

/**** Measurements procedures according to INTEL white paper

	"How to benchmark code execution times on INTEL IA-32 and IA-64"

*****/

void quicksort(uint_fast64_t* t, int n)
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

uint_fast64_t *quartiles(uint_fast64_t *tab, uint_fast64_t size)
{
	uint_fast64_t *result ;
	uint_fast64_t aux ;
	
	result = malloc(3*sizeof(uint_fast64_t));
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

static inline uint_fast64_t cpucyclesStart(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n    "
			"RDTSC\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint_fast64_t)lo)^(((uint_fast64_t)hi)<<32);
}

static inline uint_fast64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			"CPUID\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint_fast64_t)lo)^(((uint_fast64_t)hi)<<32);
}

static int_fast64_t randomint64(void)
{
	return (((int_fast64_t)rand() ^ rand()) << 32) | ((int_fast64_t)rand() ^ rand());
}

uint_fast8_t translatequaltest(poly A, poly AA)
{
	__int128 tmpp;
	uint_fast64_t Tlo[N], aux[N];
	
	for(register uint_fast16_t i = 0; i < N; i++)
		Tlo[i] = translaTlo[i] + (A->t[i] - AA->t[i]);
	
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		aux[i] = (uint_fast64_t)B1[0][i]*Tlo[0];
		for (int j = 1; j < N; j++)
			aux[i] += (uint_fast64_t)B1[j][i]*Tlo[j];
	}
	
	int_fast64_t chk = 0;
	register uint_fast16_t i = 0;
	while((i < N) && chk == 0)
	{
		tmpp = (__int128)B[0][i]*aux[0];
		for (register uint_fast16_t j = 1; j < N; j++)
			tmpp += (__int128)B[j][i]*aux[j];
		
		chk = translaThi[i] + (tmpp >> 64) + 1;
		i += 1;
	}
	return chk == 0;
}

uint_fast8_t floatbabaiequaltest(poly A, poly AA)
{
	int_fast64_t V[N];
	int_fast8_t T[N];
	float temporum;
	for(register uint_fast16_t i = 0; i < N; i++)
		V[i] = A->t[i] - AA->t[i];
	
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		temporum = floL1[0][i]*V[0];
		for (register uint_fast16_t j = 1; j < N; j++)
			temporum += floL1[j][i]*V[j];
		T[i] = temporum + 0.5*(1 - 2*(temporum < 0));
	}
	
	register uint_fast16_t i = 0;
	register uint_fast8_t chk = 0;
	//for(register uint_fast16_t i = 0; i < N; i++)
	while((i < N) && chk == 0)
	{
		//S = B[0][i]*T[0];
		for (register uint_fast16_t j = 0; j < N; j++)
			V[i] -= B[j][i]*T[j];
		chk = (V[i]);
		i += 1;
	}
	return !chk;
}

/*uint_fast8_t avxfloatbabaiequaltest(poly A, poly AA)
{
	int_fast64_t V[N];
	int_fast8_t T[N];
	//float temporum;
	fvector temporum, vvect = {0}, sumvect1, sumvect2;
	float temporium;
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		V[i] = A->t[i] - AA->t[i];
		vvect.f[i] = V[i];
	}
	
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		//temporum = floL1[0][i]*V[0];
		//for (register uint_fast16_t j = 1; j < N; j++)
		//	temporum += floL1[j][i]*V[j];
		//T[i] = temporum + 0.5*(1 - 2*(temporum < 0));
		temporum.m = _mm256_mul_ps(vvect.m, fvectL1[i].m);
		temporium = temporum.f[0];
		for(register uint_fast16_t j = 1; j < N; j++)
			temporium += temporum.f[j];
		//sumvect1.m = _mm256_hadd_ps(temporum.m, temporum.m);
		//sumvect2.m = _mm256_hadd_ps(sumvect1.m, sumvect1.m);
		//temporium = sumvect1.f[0] + sumvect1.f[1] + temporum.f[5];
		T[i] = temporium + 0.5*(1 - 2*(temporium < 0));
	}
	
	register uint_fast16_t i = 0;
	register uint_fast8_t chk = 0;
	//for(register uint_fast16_t i = 0; i < N; i++)
	while((i < N) && chk == 0)
	{
		//S = B[0][i]*T[0];
		for (register uint_fast16_t j = 0; j < N; j++)
			V[i] -= B[j][i]*T[j];
		chk = (V[i]);
		i += 1;
	}
	return !chk;
}*/

uint_fast8_t ratbabaiequaltest(poly A, poly AA)
{
	int_fast64_t V[N], S;
	int8_t T[N];
	__int128 temporum;
	int8_t *movt = &temporum;
	for(register uint_fast16_t i = 0; i < N; i++)
		V[i] = A->t[i] - AA->t[i];
	
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		temporum = (__int128)ratL1[0][i]*V[0];
		for (register uint_fast16_t j = 1; j < N; j++)
			temporum += (__int128)ratL1[j][i]*V[j];
		T[i] = (movt[9]+(1 - 2*(movt[9] < 0)))/4;
		//printf("%f\n", movt[9]/4.0);
	}
	
	for(register uint_fast16_t i = 0; i < N; i++)
	{
		S = B[0][i]*T[0];
		for (register uint_fast16_t j = 1; j < N; j++)
			S += B[j][i]*T[j];
		if(S != V[i])
			return 0;
	}
	return 1;
}

void inequalitybench(uint_fast64_t retcycles[3], uint_fast8_t (*equalitytest1)(poly, poly), uint_fast8_t (*equalitytest2)(poly, poly), uint_fast8_t (*equalitytest3)(poly, poly), const uint_fast64_t W)
{
	uint_fast64_t *cycles = (uint_fast64_t *)calloc(NTEST,sizeof(uint_fast64_t)), *statTimer;
	uint_fast64_t timermin, timermax, meanTimermin = 0, meanTimermax = 0, t1, t2, diff_t;
	poly a, b;
	init_polys(N, &a, &b, NULL);
	uint_fast64_t c1 = 0, c2 = 0, c3 = 0;
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(a);
		randpoly(b);
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c1 += equalitytest1(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c1 += equalitytest1(a, b);
				c1 += equalitytest1(b, a);
				c1 += equalitytest1(a, b);
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
		retcycles[0] += statTimer[1];
		free(statTimer);
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c2 += equalitytest2(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c2 += equalitytest2(a, b);
				c2 += equalitytest2(b, a);
				c2 += equalitytest2(a, b);
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
		retcycles[1] += statTimer[1];
		free(statTimer);
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c3 += equalitytest3(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c3 += equalitytest3(a, b);
				c3 += equalitytest3(b, a);
				c3 += equalitytest3(a, b);
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
		retcycles[2] += statTimer[1];
		free(statTimer);
	}
	
	
	printf("check: %ld %ld %ld\n", c1, c2, c3);
	
	free(cycles);
	free_polys(a, b, NULL);
	retcycles[0] = retcycles[0]/NSAMPLES/W; // We divide by W since we measured W calls.
	retcycles[1] = retcycles[1]/NSAMPLES/W;
	retcycles[2] = retcycles[2]/NSAMPLES/W;
}

void equalitybench(uint_fast64_t retcycles[3], uint_fast8_t (*equalitytest1)(poly, poly), uint_fast8_t (*equalitytest2)(poly, poly), uint_fast8_t (*equalitytest3)(poly, poly), const uint_fast64_t W)
{
	uint_fast64_t *cycles = (uint_fast64_t *)calloc(NTEST,sizeof(uint_fast64_t)), *statTimer;
	uint_fast64_t timermin, timermax, meanTimermin = 0, meanTimermax = 0, t1, t2, diff_t;
	poly a, b;
	init_polys(N, &a, &b, NULL);
	uint_fast64_t c1 = 0, c2 = 0, c3 = 0;
	int_fast64_t spash = 0;
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly(a);
		
		for(int i = 0; i < N; i++)
			b->t[i] = a->t[i];
		
		for(int i = 0; i < N; i++)
		{
			spash = (randomint64() % (2*N)) - N + 1;
			if(spash)
				for(int j = 0; j < N; j++)
					b->t[j] += B[i][j] * spash;
		}
		
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c1 += equalitytest1(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c1 += equalitytest1(a, b);
				c1 += equalitytest1(a, a);
				c1 += equalitytest1(b, b);
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
		retcycles[0] += statTimer[1];
		free(statTimer);
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c2 += equalitytest2(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c2 += equalitytest2(a, b);
				c2 += equalitytest2(a, a);
				c2 += equalitytest2(b, b);
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
		retcycles[1] += statTimer[1];
		free(statTimer);
		timermin = (uint_fast64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint_fast64_t));
		for(int i=0;i<NTEST;i++)
		{
		// Here we "heat" the cache memory.
			c3 += equalitytest3(a, b);
		}
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(uint_fast64_t soak=0; soak < W/3; soak++)
			{
				c3 += equalitytest3(a, b);
				c3 += equalitytest3(a, a);
				c3 += equalitytest3(b, b);
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
		retcycles[2] += statTimer[1];
		free(statTimer);
	}
	
	
	printf("check: %ld %ld %ld\n", c1, c2, c3);
	
	free(cycles);
	free_polys(a, b, NULL);
	retcycles[0] = retcycles[0]/NSAMPLES/W; // We divide by W since we measured W calls.
	retcycles[1] = retcycles[1]/NSAMPLES/W;
	retcycles[2] = retcycles[2]/NSAMPLES/W;
}

int main()
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	/*int_fast64_t tmp;
	sixtyfourshift(&tmp, 137);
	printf("%ld\n", tmp);
	thirtytwoshift(&tmp, 137);
	printf("%ld\n", tmp);
	floamul(&tmp, 137);
	printf("%ld\n", tmp);*/
	//printf("%ld\n%ld\n%ld\n", do_approxbench(sixtyfourshift, 3003), do_approxbench(thirtytwoshift, 3003), do_approxbench(floamul, 3003));
	
	poly A, AA;
	init_polys(N, &A, &AA, NULL);
	
	randpoly(A);
	
	poly_print(A);
	
	randpoly(AA);
	
	for(int i = 0; i < N; i++)
		AA->t[i] = A->t[i];
	
	int_fast64_t spash; // = randomint64() & ((1<<N) - 1);
	
	for(int i = 0; i < N; i++)
	{
		spash = (randomint64() % (2*N)) - N + 1;
		//if ((1<<i) & spash)
		if(spash)
			for(int j = 0; j < N; j++)
				AA->t[j] += B[i][j] * spash;
	}
	
	poly_print(AA);
	
	printf("%d\n", translatequaltest(A, AA));
	printf("%d\n", floatbabaiequaltest(A, AA));
	printf("%d\n", ratbabaiequaltest(A, AA));
	
	
	/*_poly display;
	display.deg = N;
	display.t = V;
	poly_print(&display);
	display.t = T;
	poly_print(&display);
	display.t = S;
	poly_print(&display);*/
	uint_fast64_t cycles[3] = {0}, ecycles[3] = {0};
	const uint_fast64_t rep = 99;
	printf("equality:\n");
	equalitybench(ecycles, floatbabaiequaltest, ratbabaiequaltest, translatequaltest, rep);
	printf("| babai float | %ld |\n| babai int   | %ld |\n| translation | %ld |\n", ecycles[0], ecycles[1], ecycles[2]);
	printf("inequality:\n");
	inequalitybench(cycles, floatbabaiequaltest, ratbabaiequaltest, translatequaltest, rep);
	printf("| babai float | %ld |\n| babai int   | %ld |\n| translation | %ld |\n", cycles[0], cycles[1], cycles[2]);
	
	free_polys(A, AA, NULL);
	return 0;
}

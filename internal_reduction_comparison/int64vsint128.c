#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
//#include <immintrin.h>

#define NTEST 501
#define NSAMPLES 1001

#define N 9

uint8_t indices[N*N] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 0, 2, 3, 4, 5, 6, 7, 8, 0, 1, 3, 4, 5, 6, 7, 8, 0, 1, 2, 4, 5, 6, 7, 8, 0, 1, 2, 3, 5, 6, 7, 8, 0, 1, 2, 3, 4, 6, 7, 8, 0, 1, 2, 3, 4, 5, 7, 8, 0, 1, 2, 3, 4, 5, 6, 8, 0, 1, 2, 3, 4, 5, 6, 7};

// NTEST*NSAMPLES must be odd
// it's easier to compute median value

#ifdef RDPMC_ALLOWED

unsigned long rdpmc_instructions(void) { return 1;}

#else

inline static unsigned long rdpmc_instructions(void)
{
	unsigned a, d, c;
	
	c = (1<<30);
	__asm__ __volatile__("rdpmc" : "=a" (a), "=d" (d) : "c" (c));
	
	return ((unsigned long)a) | (((unsigned long)d) << 32);;
}

#endif

enum intmultmode{S6464, S64128, S128128, S646464, S3232, S323232, S1616};

void __print128(register const __int128 Val)
{
	int64_t hi = Val >> 64;
	uint64_t lo = Val;
	printf("%lx%016lx", hi, lo);
}

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

__int128 randomint128(void)
{
	return ((__int128)randomint64() << 64) | randomint64();
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

void intmult6464(__int128* res, const int64_t* op1, const int64_t* op2)
{
	//*res += (__int128) op1 * op2;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += (__int128) op1[i] * op2[j];
}
void intmult64128(__int128* res, const int64_t* op1, const __int128* op2)
{
	//*res += (__int128) op2 * op1;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += (__int128) op2[i] * op1[j];
}
void intmult128128(__int128* res, const __int128* op1, const __int128* op2)
{
	//*res += (__int128) op1 * op2;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += (__int128) op1[i] * op2[j];
}
void intmult646464(int64_t* res, const int64_t* op1, const int64_t* op2)
{
	//*res += op1 * op2;
	//long long unsigned int dump;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] +=  op1[i] * op2[j];
			//res[indices[i+j]] += _mulx_u64(op1[i], op2[i], &dump);
}
void intmult3232(int64_t* res, const int32_t* op1, const int32_t* op2)
{
	//*res = (int64_t)op1 * op2;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += (int64_t) op1[i] * op2[j];
}
void intmult323232(int32_t* res, const int32_t* op1, const int32_t* op2)
{
	//*res = op1 * op2;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += op1[i] * op2[j];
}
void intmult1616(int32_t* res, const int16_t* op1, const int16_t* op2)
{
	//*res = (int32_t)op1 * op2;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			res[indices[i+j]] += (int32_t) op1[i] * op2[j];
}

void do_bench(uint64_t retcycles[3], const uint64_t W, const enum intmultmode MODE)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin, timermax, meanTimermin =0, medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	/*__int128 a = 0, b = 0, c = randomint128();
	int32_t ku = 0, kul = 0, kan = 0;
	int64_t Gucumatz = 0;
	int16_t quetzal = 0, coatl = 0;*/
	__int128 a[N] = {0}, b[N] = {0}, c[N];
	for(int i = 0; i < N; i++)
		c[i] = randomint128();
	int32_t ku[N] = {0}, kul[N] = {0}, kan[N] = {0};
	int64_t Gucumatz[N] = {0};
	int16_t quetzal[N] = {0}, coatl[N] = {0};
	__int128 tmptoprint128 = 0; int64_t tmptoprint64 = 0; int32_t tmptoprint32 = 0;
	for(int i = 0; i < NTEST; i++)
	{
	// Here we "heat" the cache memory.
		switch(MODE)
		{
			case S6464:
				for(int j = 0; j < N; j++)
				{
					a[j] = randomint128();
					b[j] = randomint128();
				}
				intmult6464(c, (int64_t*)a, (int64_t*)b);
				break;
			case S64128:
				for(int j = 0; j < N; j++)
				{
					a[j] = randomint128();
					b[j] = randomint128();
				}
				intmult64128(c, (int64_t*)a, b);
				break;
			case S128128:
				for(int j = 0; j < N; j++)
				{
					a[j] = randomint128();
					b[j] = randomint128();
				}
				intmult128128(c, a, b);
				break;
			case S646464:
				for(int j = 0; j < N; j++)
				{
					a[j] = randomint128();
					b[j] = randomint128();
				}
				intmult646464((int64_t*)c, (int64_t*)a, (int64_t*)b);
				break;
			case S3232:
				for(int j = 0; j < N; j++)
				{
					ku[j] = randomint64();
					kul[j] = randomint64();
				}
				intmult3232(Gucumatz, ku, kul);
				break;
			case S323232:
				for(int j = 0; j < N; j++)
				{
					ku[j] = randomint64();
					kul[j] = randomint64();
				}
				intmult323232(kan, ku, kul);
				break;
			case S1616:
				for(int j = 0; j < N; j++)
				{
					quetzal[j] = randomint64();
					coatl[j] = randomint64();
				}
				intmult1616(kan, quetzal, coatl);
				break;
			default:
				exit(1);
		}
		for(int j = 0; j < N; j++)
			tmptoprint128 += (__int128)c[j];
		for(int j = 0; j < N; j++)
			tmptoprint64 += (int64_t)Gucumatz[j];
		for(int j = 0; j < N; j++)
			tmptoprint32 += (int32_t)kan[j];
	}
	
	for(int i = 0; i < NSAMPLES; i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		switch(MODE)
		{
			case S3232:
			case S323232:
				for(int j = 0; j < N; j++)
				{
					ku[j] = randomint64();
					kul[j] = randomint64();
				}
				break;
			case S1616:
				for(int j = 0; j < N; j++)
				{
					quetzal[j] = randomint64();
					coatl[j] = randomint64();
				}
				break;
			default:
				for(int j = 0; j < N; j++)
				{
					a[j] = randomint128();
					b[j] = randomint128();
				}
		}
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			printf("\b%d\t%d\r", i, j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(unsigned int soak = 0; soak < W/2; soak++)
			{
				switch(MODE)
				{
					case S6464:
						intmult6464(c, (int64_t*)a, (int64_t*)b);
						break;
					case S64128:
						intmult64128(c, (int64_t*)a, b);
						break;
					case S128128:
						intmult128128(c, a, b);
						break;
					case S646464:
						intmult646464((int64_t*)c, (int64_t*)a, (int64_t*)b);
						break;
					case S3232:
						intmult3232(Gucumatz, ku, kul);
						break;
					case S323232:
						intmult323232(kan, ku, kul);
						break;
					case S1616:
						intmult1616(kan, quetzal, coatl);
						break;
					default:
						exit(1);
				}
				switch(MODE)
				{
					case S6464:
					case S64128:
					case S128128:
					case S646464:
						for(int j = 0; j < N; j++)
							b[j] = c[j];
						break;
					case S3232:
						for(int j = 0; j < N; j++)
							kul[j] = Gucumatz[j];
						break;
					case S323232:
						for(int j = 0; j < N; j++)
							kul[j] = kan[j];
						break;
					case S1616:
						for(int j = 0; j < N; j++)
							coatl[j] = kan[j];
						break;
					default:
						exit(1);
				}
				switch(MODE)
				{
					case S6464:
						intmult6464(c, (int64_t*)a, (int64_t*)b);
						break;
					case S64128:
						intmult64128(c, (int64_t*)a, b);
						break;
					case S128128:
						intmult128128(c, a, b);
						break;
					case S646464:
						intmult646464((int64_t*)c, (int64_t*)a, (int64_t*)b);
						break;
					case S3232:
						intmult3232(Gucumatz, ku, kul);
						break;
					case S323232:
						intmult323232(kan, ku, kul);
						break;
					case S1616:
						intmult1616(kan, quetzal, coatl);
						break;
					default:
						exit(1);
				}
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
			for(int j = 0; j < N; j++)
				tmptoprint128 += (__int128)c[j];
			for(int j = 0; j < N; j++)
				tmptoprint64 += (int64_t)Gucumatz[j];
			for(int j = 0; j < N; j++)
				tmptoprint32 += (int32_t)kan[j];
		}
		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles,NTEST);
		medianTimer += statTimer[1];
		free(statTimer);
	}
	
	switch(MODE)
	{
		case S6464:
		case S64128:
		case S128128:
		case S646464:
			printf("0x"); __print128(tmptoprint128); printf("\r");
			break;
		case S3232:
			printf("0x%lx\r", tmptoprint64);
			break;
		case S1616:
		case S323232:
			printf("0x%x\r", tmptoprint32);
			break;
		default:
			exit(1);
	}
	
	free(cycles);
	// We divide by W since we measured W calls.
	printf("                                          \r");
	retcycles[0] = meanTimermin/NSAMPLES;
	retcycles[1] = meanTimermax/NSAMPLES;
	retcycles[2] = medianTimer/NSAMPLES;
}


int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles[3];
	
	const uint64_t W = 1000;
	
	do_bench(cycles, W, S6464);
	//printf("S6464 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC128 = A64 x B64 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S64128);
	//printf("S64128 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC128 = A64 x B128 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S128128);
	//printf("S128128 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC128 = A128 x B128 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S646464);
	//printf("S646464 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC64 = A64 x B64 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S3232);
	//printf("S3232 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC64 = A32 x B32 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S323232);
	//printf("S323232 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC32 = A32 x B32 %.3lf\n", (double)cycles[2]/W/N/N);
	
	do_bench(cycles, W, S1616);
	//printf("S1616 (%ld, %ld, %ld)\n", cycles[0], cycles[1], cycles[2]);
	printf("\bC32 = A16 x B16 %.3lf\n", (double)cycles[2]/W/N/N);
	
	return 0;
}


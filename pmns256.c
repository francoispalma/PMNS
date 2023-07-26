#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


#define NTEST 501
#define NSAMPLES 1001


#include "pmns256.h"

#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 34")

#define add_overflow(X, Y) __builtin_add_overflow(*(X), Y, X)

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

void randpoly256(restrict poly256 P)
{
	for(int i = 0; i < N; i++)
	{
		P->lo[i] = randomint64();
		P->midlo[i] = randomint64();
		P->midhi[i] = randomint64();
		P->hi[i] = randomint64() & ((1ULL << (RHO - 192)) - 1);
	}
}

void pmns256_montg_mult(restrict poly256 res, const restrict poly256 A,
	const restrict poly256 B)
{
	// Function that multiplies A by B using the Montgomery CIOS method in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will be in the PMNS also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	int64_t tmplo[N];
	uint64_t Shi[N];
	unsigned __int128 tmp, tmp2, tmp3, aux, aux2, aux3;
	unsigned __int128 Slo[N], Tlo[N], Tmidlo[N], Tmidhi[N], Thi[N];
	
	// Shi, Slo = Alo * Blo
	
	for(int i = 0; i < N - 1; i++)
	{
		Slo[i] = (unsigned __int128) A->lo[i + 1] * B->lo[N - 1];
		Shi[i] = 0;
		for(int j = 2; j < N - i; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[i + j] * B->lo[N - j]);
		}
		aux = (unsigned __int128) LOW(Slo[i]) * LAMBDA;
		tmp = (unsigned __int128) HI(Slo[i]) * LAMBDA + HI(aux);
		Slo[i] = ((unsigned __int128) tmp << 64) | LOW(aux);
		Shi[i] = (unsigned __int128) Shi[i] * LAMBDA + HI(tmp);
	}
	Slo[N - 1] = 0;
	Shi[N - 1] = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[j] * B->lo[i - j]);
		}
	}
	
	// tmplo = Slo * M1Lo
	
	for(int i = 0; i < N; i++)
	{
		tmplo[i] = Slo[0] * B1lo[0][i];
		
		for(int j = 1; j < N; j++)
			tmplo[i] += Slo[j] * B1lo[j][i];
	}
	
	// Thi, Tmidhi, Tmidlo, Tlo = M * tmplo
	
	for(int i = 0; i < N; i++)
	{
		tmp = (unsigned __int128) Blo[0][i] * tmplo[0];
		tmp2 = (unsigned __int128) Bmidlo[0][i] * tmplo[0];
		tmp3 = (unsigned __int128) Bmidhi[0][i] * tmplo[0];
		Thi[i] = (__int128) Bhi[0][i] * tmplo[0] + HIGH(tmp3);
		Tmidhi[i] = LOW(tmp3) + HIGH(tmp2);
		Tmidlo[i] = LOW(tmp2) + HIGH(tmp);
		Tlo[i] = LOW(tmp);
		for(int j = 1; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[j][i] * tmplo[j];
			tmp2 = (unsigned __int128) Bmidlo[j][i] * tmplo[j];
			tmp3 = (unsigned __int128) Bmidhi[j][i] * tmplo[j];
			Thi[i] += (__int128) Bhi[j][i] * tmplo[j] + HIGH(tmp3);
			Tmidhi[i] += LOW(tmp3) + HIGH(tmp2);
			Tmidlo[i] += LOW(tmp2) + HIGH(tmp);
			Tlo[i] += LOW(tmp);
		}
	}
	
	// Thi, Tmidhi, Tmidlo, Tlo += Blo * A
	// We already have Alo * Blo in Slo and Shi
	// In the process we also get Tmidlo += Bmidlo * Alo
	
	for(int i = 0; i < N - 1; i++)
	{
		aux = (__int128) A->midlo[i + 1] * B->lo[N - 1];
		aux2 = (__int128) A->lo[i + 1] * B->midlo[N - 1];
		aux3 = (__int128) A->midhi[i + 1] * B->lo[N - 1];
		tmp3 = (__int128) A->hi[i + 1] * B->lo[N - 1] + HIGH(aux3);
		tmp2 = LOW(aux3) + HIGH(aux) + HIGH(aux2);
		tmp = LOW(aux) + LOW(aux2);
		for(int j = 2; j < N - i; j++)
		{
			aux = (__int128) A->midlo[i + j] * B->lo[N - j];
			aux2 = (__int128) A->lo[i + j] * B->midlo[N - j];
			aux3 = (__int128) A->midhi[i + j] * B->lo[N - j];
			tmp3 += (__int128) A->hi[i + j] * B->lo[N - j] + HIGH(aux3);
			tmp2 += LOW(aux3) + HIGH(aux) + HIGH(aux2);
			tmp += LOW(aux) + LOW(aux2);
		}
		
		Tmidlo[i] += LAMBDA * tmp;
		Tmidhi[i] += Shi[i] + LAMBDA * tmp2;
		Thi[i] += LAMBDA * tmp3;
	}
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			aux = (__int128) A->midlo[j] * B->lo[i - j];
			aux2 = (__int128) A->lo[j] * B->midlo[i - j];
			aux3 = (__int128) A->midhi[j] * B->lo[i - j];
			Thi[i] += (__int128) A->hi[j] * B->lo[i - j] + HIGH(aux3);
			Tmidhi[i] += LOW(aux3) + HIGH(aux) + HIGH(aux2);
			Tmidlo[i] += LOW(aux) + LOW(aux2);
		}
		Tmidlo[i] += HI(Slo[i]) + HIGH(Tlo[i]) + 1;
	}
	
	// tmplo computes (((S % r + A_0*(B//(r**i) % r) % E) * M_prime_0 % E) % r)
	
	for(int i = 0; i < N; i++)
	{
		tmplo[i] = Tmidlo[0] * B1lo[0][i];
		
		for(int j = 1; j < N; j++)
			tmplo[i] += Tmidlo[j] * B1lo[j][i];
	}
	
	// Tlo, Thi, Tmidhi, Tmidlo += M * tmplo
	
	for(int i = 0; i < N; i++)
	{
		tmp = (unsigned __int128) Blo[0][i] * tmplo[0];
		tmp2 = (unsigned __int128) Bmidlo[0][i] * tmplo[0];
		tmp3 = (unsigned __int128) Bmidhi[0][i] * tmplo[0];
		Tlo[i] = (__int128) Bhi[0][i] * tmplo[0] + HIGH(tmp3);
		Thi[i] += LOW(tmp3) + HIGH(tmp2);
		Tmidhi[i] += LOW(tmp2) + HIGH(tmp);
		Tmidlo[i] += LOW(tmp);
		for(int j = 1; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[j][i] * tmplo[j];
			tmp2 = (unsigned __int128) Bmidlo[j][i] * tmplo[j];
			tmp3 = (unsigned __int128) Bmidhi[j][i] * tmplo[j];
			Tlo[i] += (__int128) Bhi[j][i] * tmplo[j] + HIGH(tmp3);
			Thi[i] += LOW(tmp3) + HIGH(tmp2);
			Tmidhi[i] += LOW(tmp2) + HIGH(tmp);
			Tmidlo[i] += LOW(tmp);
		}
	}
	
	// Tlo, Thi, Tmidhi, Tmidlo += Bmidlo * A
	// We already have Alo * Bmidlo in Tmidlo and Tmidhi
	// In the process we also get Tmidhi += Bmidhi * Alo
	
	for(int i = 0; i < N - 1; i++)
	{
		aux = (__int128) A->midlo[i + 1] * B->midlo[N - 1];
		aux2 = (__int128) A->lo[i + 1] * B->midhi[N - 1];
		aux3 = (__int128) A->midhi[i + 1] * B->midlo[N - 1];
		tmp3 = (__int128) A->hi[i + 1] * B->midlo[N - 1] + HIGH(aux3);
		tmp2 = LOW(aux3) + HIGH(aux) + HIGH(aux2);
		tmp = LOW(aux) + LOW(aux2);
		for(int j = 2; j < N - i; j++)
		{
			aux = (__int128) A->midlo[i + j] * B->midlo[N - j];
			aux2 = (__int128) A->lo[i + j] * B->midhi[N - j];
			aux3 = (__int128) A->midhi[i + j] * B->midlo[N - j];
			tmp3 += (__int128) A->hi[i + j] * B->midlo[N - j] + HIGH(aux3);
			tmp2 += LOW(aux3) + HIGH(aux) + HIGH(aux2);
			tmp += LOW(aux) + LOW(aux2);
		}
		
		Tmidhi[i] += LAMBDA * tmp;
		Thi[i] += LAMBDA * tmp2;
		Tlo[i] += LAMBDA * tmp3;
	}
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			aux = (__int128) A->midlo[j] * B->midlo[i - j];
			aux2 = (__int128) A->lo[j] * B->midhi[i - j];
			aux3 = (__int128) A->midhi[j] * B->midlo[i - j];
			Tlo[i] += (__int128) A->hi[j] * B->midlo[i - j] + HIGH(aux3);
			Thi[i] += LOW(aux3) + HIGH(aux) + HIGH(aux2);
			Tmidhi[i] += LOW(aux) + LOW(aux2);
		}
		Tmidhi[i] += HIGH(Tmidlo[i]) + 1;
	}
	
	// tmplo cancels the 64 next bits
	
	for(int i = 0; i < N; i++)
	{
		tmplo[i] = Tmidhi[0] * B1lo[0][i];
		
		for(int j = 1; j < N; j++)
			tmplo[i] += Tmidhi[j] * B1lo[j][i];
	}
	
	// Tmidlo, Tlo, Thi, Tmidhi += M * tmplo
	
	for(int i = 0; i < N; i++)
	{
		tmp = (unsigned __int128) Blo[0][i] * tmplo[0];
		tmp2 = (unsigned __int128) Bmidlo[0][i] * tmplo[0];
		tmp3 = (unsigned __int128) Bmidhi[0][i] * tmplo[0];
		Tmidlo[i] = (__int128) Bhi[0][i] * tmplo[0] + HIGH(tmp3);
		Tlo[i] += LOW(tmp3) + HIGH(tmp2);
		Thi[i] += LOW(tmp2) + HIGH(tmp);
		Tmidhi[i] += LOW(tmp);
		for(int j = 1; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[j][i] * tmplo[j];
			tmp2 = (unsigned __int128) Bmidlo[j][i] * tmplo[j];
			tmp3 = (unsigned __int128) Bmidhi[j][i] * tmplo[j];
			Tmidlo[i] += (__int128) Bhi[j][i] * tmplo[j] + HIGH(tmp3);
			Tlo[i] += LOW(tmp3) + HIGH(tmp2);
			Thi[i] += LOW(tmp2) + HIGH(tmp);
			Tmidhi[i] += LOW(tmp);
		}
	}
	
	// Tmidlo, Tlo, Thi, Tmidhi += Bmidhi * A
	// We already have Alo * Bmidlo in Tmidlo and Tmidhi
	// In the process we also get Thi += Bhi * Alo
	
	for(int i = 0; i < N - 1; i++)
	{
		aux = (__int128) A->midlo[i + 1] * B->midhi[N - 1];
		tmp = (__int128) A->lo[i + 1] * B->hi[N - 1] + LOW(aux);
		aux3 = (__int128) A->midhi[i + 1] * B->midhi[N - 1];
		tmp3 = (__int128) A->hi[i + 1] * B->midhi[N - 1] + HIGH(aux3);
		tmp2 = LOW(aux3) + HIGH(aux);
		for(int j = 2; j < N - i; j++)
		{
			aux = (__int128) A->midlo[i + j] * B->midhi[N - j];
			tmp += (__int128) A->lo[i + j] * B->hi[N - j] + LOW(aux);
			aux3 = (__int128) A->midhi[i + j] * B->midhi[N - j];
			tmp3 += (__int128) A->hi[i + j] * B->midhi[N - j] + HIGH(aux3);
			tmp2 += LOW(aux3) + HIGH(aux);
		}
		
		Thi[i] += LAMBDA * tmp;
		Tlo[i] += LAMBDA * tmp2;
		Tmidlo[i] += LAMBDA * tmp3;
	}
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			aux = (__int128) A->midlo[j] * B->midhi[i - j];
			Thi[i] += (__int128) A->lo[j] * B->hi[i - j] + LOW(aux);
			aux3 = (__int128) A->midhi[j] * B->midhi[i - j];
			Tmidlo[i] += (__int128) A->hi[j] * B->midhi[i - j] + HIGH(aux3);
			Tlo[i] += LOW(aux3) + HIGH(aux);
		}
		Thi[i] += HIGH(Tmidhi[i]) + 1;
	}
	
	// Last loop iteration
	
	for(int i = 0; i < N; i++)
	{
		tmplo[i] = Thi[0] * B1lo[0][i];
		
		for(int j = 1; j < N; j++)
			tmplo[i] += Thi[j] * B1lo[j][i];
	}
	
	for(int i = 0; i < N; i++)
	{
		tmp = (unsigned __int128) Blo[0][i] * tmplo[0];
		tmp2 = (unsigned __int128) Bmidlo[0][i] * tmplo[0];
		tmp3 = (unsigned __int128) Bmidhi[0][i] * tmplo[0];
		Tmidhi[i] = (__int128) Bhi[0][i] * tmplo[0] + HIGH(tmp3);
		Tmidlo[i] += LOW(tmp3) + HIGH(tmp2);
		Tlo[i] += LOW(tmp2) + HIGH(tmp);
		Thi[i] += LOW(tmp);
		for(int j = 1; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[j][i] * tmplo[j];
			tmp2 = (unsigned __int128) Bmidlo[j][i] * tmplo[j];
			tmp3 = (unsigned __int128) Bmidhi[j][i] * tmplo[j];
			Tmidhi[i] += (__int128) Bhi[j][i] * tmplo[j] + HIGH(tmp3);
			Tmidlo[i] += LOW(tmp3) + HIGH(tmp2);
			Tlo[i] += LOW(tmp2) + HIGH(tmp);
			Thi[i] += LOW(tmp);
		}
	}
	
	for(int i = 0; i < N - 1; i++)
	{
		aux = (__int128) A->midlo[i + 1] * B->hi[N - 1];
		aux2 = (__int128) A->midhi[i + 1] * B->hi[N - 1];
		aux3 = (__int128) A->hi[i + 1] * B->hi[N - 1];
		for(int j = 2; j < N - i; j++)
		{
			aux += (__int128) A->midlo[i + j] * B->hi[N - j];
			aux2 += (__int128) A->midhi[i + j] * B->hi[N - j];
			aux3 += (__int128) A->hi[i + j] * B->hi[N - j];
		}
		
		Tlo[i] += LAMBDA * aux;
		Tmidlo[i] += LAMBDA * aux2;
		Tmidhi[i] += LAMBDA * aux3;
	}
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < i + 1; j++)
		{
			Tlo[i] += (__int128) A->midlo[j] * B->hi[i - j];
			Tmidlo[i] += (__int128) A->midhi[j] * B->hi[i - j];
			Tmidhi[i] += (__int128) A->hi[j] * B->hi[i - j];
		}
		Tlo[i] += HIGH(Thi[i]) + 1;
		Tmidlo[i] += HIGH(Tlo[i]);
		Tmidhi[i] += HIGH(Tmidlo[i]);
	}
	
	for(int i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Tlo[i]);
		res->midlo[i] = LOW(Tmidlo[i]);
		res->midhi[i] = LOW(Tmidhi[i]);
		res->hi[i] = HIGH(Tmidhi[i]);
	}
}


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
	__asm__ __volatile__ (	"CPUID\n\t"
			"RDTSC\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

static inline uint64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			"CPUID\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t do_bench(void (*pmns_mult)(restrict poly256, const restrict poly256, const restrict poly256), const uint64_t W)
{
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,
	meanTimermax = 0, t1,t2, diff_t;
	stpoly256 a, b, c;
	
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		randpoly256(&a);
		randpoly256(&b);
		pmns_mult(&a, &b, &c);
	}
	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		randpoly256(&a);
		randpoly256(&b);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			// We call the function 10 times to get an accurate measurement.
			for(uint soak = 0; soak < W; soak++)
				pmns_mult(&a, &b, &c);
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
	return medianTimer/NSAMPLES/W;
}


int main()
{
	printf("%ld\n", do_bench(pmns256_montg_mult, 1));
	return 0;
}



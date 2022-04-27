#include <stdio.h>
#include <time.h>

#include "pmns128.h"
#include "utilitymp.h"

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))


static inline void mult128(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi, const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	unsigned __int128 aux, tmp;
	*Rhi = (__int128) LOW(Ahi) * LOW(Bhi);
	*Rlo = (__int128) Alo * Blo;
	aux = (__int128) Alo * LOW(Bhi);
	tmp = (__int128) *Rlo + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo > tmp);
	aux = (__int128) LOW(Ahi) * Blo;
	*Rlo = (__int128) tmp + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo < tmp);
}

static inline void umult128(unsigned __int128* Rhi, unsigned __int128* Rlo, const uint64_t Ahi, const uint64_t Alo, const uint64_t Bhi, const uint64_t Blo)
{
	unsigned __int128 aux, tmp;
	*Rhi = (__int128) LOW(Ahi) * LOW(Bhi);
	*Rlo = (__int128) Alo * Blo;
	aux = (__int128) Alo * LOW(Bhi);
	tmp = (__int128) *Rlo + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo > tmp);
	aux = (__int128) LOW(Ahi) * Blo;
	*Rlo = (__int128) tmp + (aux << 64);
	*Rhi = (__int128) *Rhi + HI(aux) + (*Rlo < tmp);
}

static inline void umultadd128(unsigned __int128* Rhi, unsigned __int128* Rlo,
	const uint64_t Ahi, const uint64_t Alo, const uint64_t Bhi,
	const uint64_t Blo)
{
	unsigned __int128 auxlo, auxhi;
	
	umult128(&auxhi, &auxlo, Ahi, Alo, Bhi, Blo);
	
	*Rlo = (__int128) *Rlo + auxlo;
	*Rhi = (__int128) *Rhi + auxhi + (*Rlo < auxlo);
}

static inline void multadd128x(__int128* Rhi, unsigned __int128* Rlo,
	const int64_t Ahi, const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	unsigned __int128 auxlo;
	__int128 auxhi;
	
	mult128(&auxhi, &auxlo, Ahi, Alo, Bhi, Blo);
	
	*Rlo = (__int128) *Rlo + auxlo;
	*Rhi = (__int128) *Rhi + auxhi + (*Rlo < auxlo);
}

void multadd128(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R;
	unsigned __int128 A0B0, A1B0, A0B1, aux4, tmplo;
	__int128 A1B1, aux1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux4 = (__int128) LOW(A0B0);
	aux3 = (__int128) HIGH(aux4) + HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);
	aux1 = (__int128) HIGH(A1B1);
	
	tmplo = *Rlo;
	*Rlo += (__int128) LOW(aux4) + (aux3 << 64);
	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);
}

void multadd128k(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R using karatsuba;
	unsigned __int128 A0B0, tmplo;
	__int128 A1B1, A1B0_A0B1l, A1B0_A0B1h, aux1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A0B0 = (__int128) Alo * Blo;
	aux3 = (__int128) (Alo - Ahi) * (Blo - Bhi);
	//A1B0_A0B1 = (__int128) A0B0 + A1B1 - aux3;
/*	printf("aux3 = ");*/
/*	__print128(aux3);*/
/*	*/
/*	printf("LOW = lambda x: x % 2**64\nHIGH = lambda x: x >> 64\nHI = lambda x: (x >> 64) % 2**64");*/
/*	printf("Ahi = 0x%lx\nAlo = 0x%lx\nBhi = 0x%lx\nBlo = 0x%lx\nA0B0 = ", Ahi, Alo, Bhi, Blo);*/
/*	*/
/*	__print128(A0B0);*/
/*	printf("A1B1 = ");*/
/*	__print128(A1B1);*/
	
	A1B0_A0B1l = (__int128) LOW(A0B0) + LOW(A1B1) - LOW(aux3);
	A1B0_A0B1h = (__int128) HI(A1B0_A0B1l) + HI(A0B0) + HIGH(A1B1) - HIGH(aux3);
	
	//__print128(A1B0_A0B1);
	//A1B0_A0B1 = (__int128) ((__int128) (A1B0_A0B1h) << 64) | LOW(A1B0_A0B1l);
/*	if(tmp != A1B0_A0B1)*/
/*	{*/
/*		__print128(tmp);*/
/*		__print128(A1B0_A0B1);*/
/*		printf("HAAAAAAAAAAA\n");*/
/*	}*/
	//__print128(A1B0_A0B1);
	//exit(0);
	
/*	aux3 = (__int128) HI(A0B0) + LOW(A1B0_A0B1);*/
/*	aux2 = (__int128) HIGH(aux3) + HIGH(A1B0_A0B1) + LOW(A1B1);*/
/*	aux1 = (__int128) HIGH(A1B1);*/
	
/*	printf("A1B0_A0B1 = ");*/
/*	__print128(A1B0_A0B1);*/
/*	printf("A1B0_A0B1l = ");*/
/*	__print128(A1B0_A0B1l);*/
/*	printf("A1B0_A0B1h = ");*/
/*	__print128(A1B0_A0B1h);*/
/*	*/
/*	printf("HI(A0B0) + LOW(A1B0_A0B1) == HI(A0B0) + LOW(A1B0_A0B1l)\n");*/
/*	printf("HIGH(aux3) + LOW(A1B0_A0B1h) + LOW(A1B1) == HIGH(aux3) + HIGH(A1B0_A0B1) + LOW(A1B1)\n");*/
/*	printf("HIGH(A1B1) == HIGH(A1B1) + HIGH(A1B0_A0B1h)\n");*/
/*	*/
/*	exit(0);*/
	
	aux3 = (__int128) HI(A0B0) + LOW(A1B0_A0B1l);
	aux2 = (__int128) HIGH(aux3) + LOW(A1B0_A0B1h) + LOW(A1B1);
	aux1 = (__int128) HIGH(A1B1);
	
	tmplo = *Rlo;
	*Rlo += (__int128) LOW(A0B0) + (aux3 << 64);
	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);
}

void m_multadd128(__int128* Rhi, unsigned __int128* Rlo, const uint64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R for mult by M use;
	unsigned __int128 A0B0, tmplo;
	__int128 A1B1, A1B0_A0B1, aux1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A0B0 = (__int128) Alo * Blo;
	aux3 = (__int128) (Alo - Ahi) * (Blo - Bhi);
	A1B0_A0B1 = (__int128) A0B0 + A1B1 - aux3;
	
	aux3 = (__int128) HI(A0B0) + LOW(A1B0_A0B1);
	aux2 = (__int128) HI(aux3) + HI(A1B0_A0B1) + LOW(A1B1);
	aux1 = (__int128) HI(A1B1);
	
	tmplo = *Rlo;
	*Rlo += (__int128) LOW(A0B0) + ((__int128)(HI(A0B0) + LOW(A1B0_A0B1)) << 64);
	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);
}

void mm1_multadd128(__int128* Rhi, unsigned __int128* Rlo, const uint64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R for mult by M or M1 use (slow);
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);
	aux1 = (__int128) HIGH(A1B1);
	
	tmplo = *Rlo;
	*Rlo += (__int128) LOW(A0B0) + (aux3 << 64);
	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);
}

void m1_multadd128(unsigned __int128* Rlo, const uint64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// multiplies A and B and adds the result to R for mult by M1 use;
	unsigned __int128 A1B0, A0B1;
	__int128 aux3;
	
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	
	aux3 = (__int128) LOW(A0B1) + LOW(A1B0);
	
	*Rlo += (__int128) Alo * Blo + (aux3 << 64);
}

static inline void mns128_mod_mult_ext_red(__int128* Rhi,
	unsigned __int128* Rlo, const restrict poly128 A, const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	unsigned __int128 aux;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128k(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		aux = (unsigned __int128) LOW(Rlo[i]) * (LAMBDA);
		aux = (unsigned __int128) HI(Rlo[i]) * (LAMBDA) + HIGH(aux);
		Rlo[i] = (__int128) Rlo[i] * (LAMBDA);
		Rhi[i] = (__int128) Rhi[i] * (LAMBDA) + HIGH(aux);
		
		for(j = 0; j < i + 1; j++)
			multadd128k(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

static inline void m_mns128_mod_mult_ext_red(__int128* Rhi, 
	unsigned __int128* Rlo, const restrict poly128 A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			mm1_multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			mm1_multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				Mhi[i - j], Mlo[i - j]);
	}
}

static inline void m1_mns128_mod_mult_ext_red(unsigned __int128* Rlo,
	const restrict poly128 A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand (in pmns128 we only care about the lower 128 bits for
	// this operation).
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			m1_multadd128(Rlo + i, A->hi[i + j], A->lo[i + j],
				M1Lambdahi[N - j], M1Lambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			m1_multadd128(Rlo + i, A->hi[j], A->lo[j],
				M1hi[i - j], M1lo[i - j]);
	}
}

static inline void mns128_montg_int_red(poly128 res, __int128* Rhi,
	unsigned __int128* Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho.
	unsigned __int128 V[N], V2[N];
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = Rlo[i];
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HI(Rlo[i]);
		V2[i] = Rhi[i];
		Rhi[i] = 0;
		Rlo[i] = 0;
	}
	
	m1_mns128_mod_mult_ext_red(Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HIGH(Rlo[i]);
		Rhi[i] = 0;
		Rlo[i] = 0;
	}
	
	m_mns128_mod_mult_ext_red(Rhi, Rlo, res);
	
	for(i = 0; i < N; i++)
	{
		Rhi[i] = (__int128) V2[i] + Rhi[i] + (V[i] != 0);
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

void convert_string_to_amns128(restrict poly128 res, const char* string)
{
	uint8_t counter;
	register uint16_t i, j;
	const unsigned __int128 rho = ((__int128)1) << RHO;
	unsigned __int128 limb, Rlo[N] = {0};
	__int128 Rhi[N] = {0}, tmp[N] = {0};
	poly stok;
	init_poly(2 * N, &stok);
	
	convert_string_to_poly(&stok, string);
	
	if(stok->deg > 2 * N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	counter = 0;
	for(i = 0; i < N - 1; i++)
	{
		limb = (((unsigned __int128) stok->t[2 * i + 1]) << 64)
		 | ((uint64_t) stok->t[2 * i]);
		tmp[i] |= (limb << counter) & (rho - 1);
		tmp[i + 1] |= (limb >> (RHO - counter));
		counter = (counter + 128 - RHO) % RHO;
	}
	limb = (unsigned __int128) ((unsigned __int128) stok->t[2 * i]) |
			(((unsigned __int128) stok->t[2 * i + 1]) << 64);
	tmp[i] |= (limb << counter) & (rho - 1);
	
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			multadd128(Rhi + j, Rlo + j, HIGH(tmp[i]), LOW(tmp[i]),
				__Pihi__[i][j], __Pilo__[i][j]);
	
	mns128_montg_int_red(res, Rhi, Rlo);
	
end:
	free_poly(stok);
}

inline void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	mns128_montg_int_red(res, Rhi, Rlo);
}

static inline int64_t randomint64(void)
{
	// Function to generate a random 64 bit number.
	return (((int64_t)rand() + rand()) << 32) ^ ((int64_t)rand() + rand());
}

static inline int64_t __modrhohi(int64_t param)
{
	// Utility function to get a high part for a random poly128.
	return param & ((1ULL<<(RHO - 64)) - 1);
}

void randpoly128(poly128 P)
{
	// Generates a random poly128 with appropriate, lower than rho high part.
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->lo[i] = randomint64();
		if(RHO > 64) P->hi[i] = __modrhohi(randomint64());
	}
}

void __main__(void)
{
/*	int64_t Ahi = 0x5c096e6b558ba549, Bhi = 0x5918460d4a05af9c;*/
/*	uint64_t  Alo = 0xf9dadbd740482391, Blo = 0xc4703ac88b18e4c3;*/
/*	*/
/*	const char RES[] = */
/*	"0x2008017506164ba2bbc36636e65b8f4caca8fc6b0e76b6a3dd70ff9147383b73";*/
	
/*	int64_t Ahi = 0xfc096e6b558ba549, Bhi = 0xf918460d4a05af9c;*/
/*	uint64_t  Alo = 0xf9dadbd740482391, Blo = 0xc4703ac88b18e4c3;*/
/*	*/
/*	const char RES[] = */
/*	"0x1b5dc7ca3fcbcc34673dbafa172c2d2ca8fc6b0e76b6a3dd70ff9147383b73";*/
	
/*	int64_t Ahi = 0xffffffffffffffff, Bhi = 0xffffffffffffffff;*/
/*	uint64_t  Alo = 0xffffffffffffffff, Blo = 0xffffffffffffffff;*/
/*	*/
/*	const char RES[] = */
/*	"0x0000000000000000000000000000000000000000000000000";*/
	
/*	int64_t Ahi = 0x7e4deecd7fb21228, Bhi = -0x11fb30a7f4f269dd;*/
/*	uint64_t  Alo = 0xbe23c803cef47a1, Blo = 0xc97d66ba9ffee37;*/
/*	*/
/*	const char RES[] = */
/*	"0xf720e4b9bc91c7b14365a9592434dee8edbb60fd615d17fd17ab67202e5f1197";*/
	
	int64_t Ahi = 0xead7e6d770b694af, Bhi = 0xa3c338cf5fdb7ee7;
	uint64_t  Alo = 0xcd3ad86ad0c870de, Blo = 0x83ed23c6860ab850;
	
	const char RES[] =
		"0x23fbe09528a6b4ada3fa9544c0b0b34f51d40d5578d6c57b14af738e2a3c5388";
	
/*	__int128 R1 = 0;*/
/*	unsigned __int128 R2 = 0;*/
	
	__int128 R1 = (__int128) ((__int128) 0x8dc15207aec5564b << 64) + 0x4ba591887a57f6ed;
	unsigned __int128 R2 = (__int128) ((__int128) 0x4b5f46c669762ca7 << 64) + 0x4f4e53803dcd7e28;
	
/*	__int128 R1 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;*/
/*	unsigned __int128 R2 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;*/
	
	multadd128(&R1, &R2, Ahi, Alo, Bhi, Blo);
	
	printf("%s\n", RES);
	
	printf("0x%lx%016lx%016lx%016lx\n", HIGH(R1), LOW(R1), HIGH(R2), LOW(R2));
	
	printf("\n\n");
	
	poly128 A;
	init_poly128(N, &A);
	const char a[] = "bf6dc9f34905d4ccea18b34313d7f22412795efa0161f7ebcd5912a900ea7d255661bb894729e4fd85a477d3c575f3e97fcd1e6e2fd01d5317724f38def3c7f944162bb4ae4dcd5b1522efca1f3713a927c91f1113096ced7585edf7fef8cc9334dc56e8483a3c49f4a0fb9bb73c00b8b00e3d11435184eacbd45dd38fcbcadd";
	
	convert_string_to_amns128(A, a);
	
	p128_print(A);
	
	free_poly128(A);
}

void __benchmult__(void)
{
	uint64_t sum1, sum2, c = clock(), lo1, lo2;
	int64_t hi1, hi2;
	__int128 dummy11 = 0, dummy21 = 0;
	unsigned __int128 dummy12 = 0, dummy22 = 0;
	
	srand((unsigned) (time(&hi1)));
	
	sum1 = 0;
	sum2 = 0;
	for(int i = 0; i < 100000; i++)
	{
		hi1 = randomint64();
		hi2 = randomint64();
		lo1 = randomint64();
		lo2 = randomint64();
		c = clock();
		multadd128x(&dummy11, &dummy12, hi1, lo1, hi2, lo2);
		sum1 += clock() - c;
		c = clock();
		multadd128(&dummy21, &dummy22, hi1, lo1, hi2, lo2);
		sum2 += clock() - c;
	}
	sum1 = 0;
	sum2 = 0;
	for(int i = 0; i < 100000; i++)
	{
		hi1 = randomint64();
		hi2 = randomint64();
		lo1 = randomint64();
		lo2 = randomint64();
		c = clock();
		multadd128x(&dummy11, &dummy12, hi1, lo1, hi2, lo2);
		sum1 += clock() - c;
		c = clock();
		multadd128(&dummy21, &dummy22, hi1, lo1, hi2, lo2);
		sum2 += clock() - c;
	}
	printf("1: %ld\n2: %ld\n", sum1, sum2);
	__print128(dummy11);
	__print128(dummy21);
	__print128(dummy12);
	__print128(dummy22);
}

void __multchecks__(void)
{
	poly128 a, b, c;
	init_poly128s(N, &a, &b, &c, NULL);
	int64_t seed;
	
	srand((unsigned) (time(&seed)));
	
	for(int i = 0; i < 1000000; i++)
	{
		randpoly128(a);
		randpoly128(b);
		p128_print(a);
		p128_print(b);
		amns128_montg_mult(c, a, b);
		p128_print(c);
	}
	
	free_poly128s(a, b, c, NULL);
}

void __multbench__(void)
{
	uint64_t c1 = 0, sum;
	poly128 a, b, c, soak1, soak2;
	init_poly128s(N, &a, &b, &c, &soak1, &soak2, NULL);
	
	srand((unsigned) (time(((int64_t*)(&c1)))));
	
	c1 = clock();
	
	randpoly128(soak2);
	soak2->lo[0] += Gi[0].t[0];
	soak2->lo[0] += __P__.t[0];
	
	for(int i = 0; i < 1000; i++)
	{
		randpoly128(a);
		randpoly128(b);
		amns128_montg_mult(c, a, b);
		amns128_montg_mult(soak1, c, soak2);
		amns128_montg_mult(soak2, c, soak1);
	}
	
	sum = 0;
	for(int i = 0; i < 10000; i++)
	{
		randpoly128(a);
		randpoly128(b);
		c1 = clock();
		amns128_montg_mult(c, a, b);
		sum += clock() - c1;
		amns128_montg_mult(soak1, c, soak2);
		amns128_montg_mult(soak2, c, soak1);
	}
	
	printf("%ld\n", sum);
	
	free_poly128s(a, b, c, soak1, soak2, NULL);
}

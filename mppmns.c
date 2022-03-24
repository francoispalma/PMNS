#include <stdio.h>

#include "mppmns.h"
#include "utilitymp.h"

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((uint64_t)(X>>64))


/*
generate P
K = GF(P)
polK.<X> = K[]
P = X^n - lambda
factor(P)
lone factor => nth root of lambda.
*/

void multadd128(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	__int128 aux1, aux3;
	unsigned __int128 aux2;
	
	// karatsuba
	aux1 = (__int128) Ahi * Bhi;
	aux2 = (__int128) Alo * Blo;
	aux3 = (__int128) ((__int128) Ahi + Alo) * ((__int128) Bhi + Blo)
		- aux1 - aux2;
	
	*Rlo += aux2 + (((__int128) LOW(aux3)) << 64);
	*Rhi += (*Rlo < aux2) + aux1 + ((__int128) HIGH(aux3));
}

static inline void mns128_mod_mult_ext_red(__int128* Rhi,
	unsigned __int128* Rlo, const restrict poly128 A, const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		Rlo[i] *= LAMBDA;
		Rhi[i] *= LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

static inline void m_mns128_mod_mult_ext_red(__int128* R, const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	/*register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) ((uint64_t)A->t[i + j]) * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) ((uint64_t)A->t[j]) * M[i - j];
	}*/
}

static inline void m1_mns128_mod_mult_ext_red(__int128* R, const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand.
	
	/*register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) ((uint64_t)A->t[i + j]) * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) ((uint64_t)A->t[j]) * M1[i - j];
	}*/
}

static inline void mns128_montg_int_red(poly128 res, const __int128* R)
{
	/*uint64_t V[N], V2[N], T[N], T2[N];
	
	for(int i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	m1_mns_mod_mult_ext_red(R, res);
	
	for(int i = 0; i < N; i++)
	{
		res->t[i] = R[i];
		R[i] = 0;
	}
	
	m_mns_mod_mult_ext_red(R, res);
	
	for(int i = 0; i < N; i++)
	{
		T[i] = R[i];
		T2[i] = (R[i] >> 64);
		
		T[i] = V[i] + T[i];
		res->t[i] = V2[i] + T2[i] + (T[i] < V[i]);
	}*/
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
			multadd128(Rhi + i, Rlo + i, HIGH(tmp[i]), LOW(tmp[i]),
				__Pihi__[i][j], __Pilo__[i][j]);
	
	mns128_montg_int_red(res, Rhi);
	
end:
	free_poly(stok);
}


int main(void)
{
	int64_t Ahi = 0x5c096e6b558ba549, Bhi = 0x5918460d4a05af9c;
	uint64_t  Alo = 0xf9dadbd740482391, Blo = 0xc4703ac88b18e4c3;
	
	const char RES[] = 
	"0x2008017506164ba2bbc36636e65b8f4caca8fc6b0e76b6a3dd70ff9147383b73";
	
	__int128 R1 = 0;
	unsigned __int128 R2 = 0;
	
	multadd128(&R1, &R2, Ahi, Alo, Bhi, Blo);
	
	printf("%s\n", RES);
	
	printf("0x%lx%016lx%016lx%016lx\n", HIGH(R1), LOW(R1), HIGH(R2), LOW(R2));
	
	printf("\n\n");
	
	const char a[] = "74ff560400d0105e6381e4f7cf22ba4a3d949bbe3b03e7ec1c8aebfb02a4dedf230eef099cd1ae78adf8f142cd70ed93122a5c48c5edcba658615fa2316994dce0c84e9e54c5ae9482acdc0ed6fae84eb7e83d94016d12452ad41369e33a53a676d539439488bdc8b3462c5579a432e8b579e8af9d5b2b0b8f37856fe2de7f30";
	
	convert_string_to_amns128(NULL, a);
	
	return 0;
}

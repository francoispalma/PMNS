#include <stdio.h>

#include "mppmns.h"
#include "utilitymp.h"

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((int64_t)(X>>64))
#define HI(X) ((uint64_t)(X>>64))

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
	// multiplies A and B and adds the result to R;
	unsigned __int128 A0B0, aux4, tmplo;
	__int128 A1B1, A1B0, A0B1, aux1, aux2, aux3;
	
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux4 = (__int128) LOW(A0B0);
	aux3 = (__int128) HIGH(aux4) + HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);
	aux1 = (__int128) HI(A1B1);
	
	tmplo = *Rlo;
	*Rlo += (__int128) LOW(aux4) + (aux3 << 64);
	*Rhi += (__int128) aux2 + (aux1 << 64) + (*Rlo < tmplo);
}

void multadd128k(__int128* Rhi, unsigned __int128* Rlo, const int64_t Ahi,
	const uint64_t Alo, const int64_t Bhi, const uint64_t Blo)
{
	// karatsuba variant, non functional
	__int128 aux1, aux3, tmp = 0, carry = 0;
	unsigned __int128 aux2, tmp2;
	
	aux1 = (__int128) Ahi * Bhi;
	aux2 = (__int128) Alo * Blo;
	//aux3 = (__int128) ((__int128) Ahi + Alo) * ((__int128) Bhi + Blo) - aux2 - aux1;
	aux3 = ((__int128) Ahi + Alo);
	tmp = ((__int128) Bhi + Blo);
	carry = (aux3 >> 62) * (tmp >> 62) >> 4;
	aux3 = (__int128) aux3 * tmp;
	tmp = aux3;
	aux3 = (__int128) aux3 - aux1;
	carry -= aux3 > tmp;
	tmp = aux3;
	aux3 = (__int128) aux3 - aux2;
	carry -= aux3 > tmp;
	
	
	tmp2 = aux2 + (((__int128) LOW(aux3)) << 64);
	*Rlo += tmp2;
	*Rhi += (tmp2 < aux2) + (tmp2 > *Rlo) + aux1 + ((__int128) HIGH(aux3)) + (carry << 64);
}

static inline void mns128_mod_mult_ext_red(__int128* Rhi,
	unsigned __int128* Rlo, const restrict poly128 A, const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	unsigned __int128 aux4;
	__int128 aux3, aux2, aux1;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		
		aux4 = (__int128) LOW(Rlo[i]) * LAMBDA;
		aux3 = (__int128) HIGH(aux4) + HIGH(Rlo[i]) * LAMBDA;
		aux2 = (__int128) HIGH(aux3) + LOW(Rhi[i]) * LAMBDA;
		aux1 = (__int128) HIGH(Rhi[i]) * LAMBDA;
		
		Rlo[i] = (__int128) LOW(aux4) + (aux3 << 64);
		Rhi[i] = (__int128) aux2 + (aux1 << 64);
		
/*		Rlo[i] = (__int128) Rlo[i] * LAMBDA;*/
/*		Rhi[i] = (__int128) Rhi[i] * LAMBDA;*/
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
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
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				Mhi[i - j], Mlo[i - j]);
	}
}

static inline void m1_mns128_mod_mult_ext_red(__int128* Rhi, 
	unsigned __int128* Rlo, const restrict poly128 A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[i + j], A->lo[i + j],
				M1Lambdahi[N - j], M1Lambdalo[N - j]);
		
		printf("0x%lx%016lx%016lx%016lx, ", HIGH(Rhi[i]), LOW(Rhi[i]), HIGH(Rlo[i]), LOW(Rlo[i]));
		
		for(j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				M1hi[i - j], M1lo[i - j]);
	}
}

static inline void mns128_montg_int_red(poly128 res, __int128* Rhi,
	unsigned __int128* Rlo)
{
	unsigned __int128 V[N], V2[N], T[N], T2[N];
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = Rlo[i];
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HIGH(Rlo[i]);
		V2[i] = Rhi[i];
		Rhi[i] = 0;
		Rlo[i] = 0;
	}
	
	m1_mns128_mod_mult_ext_red(Rhi, Rlo, res);
	
	
	printf("[");
	for(i = 0; i < N; i++)
		printf("0x%lx%016lx%016lx%016lx, ", HIGH(Rhi[i]), LOW(Rhi[i]), HIGH(Rlo[i]), LOW(Rlo[i]));
	printf("]\n\n");
	
	/*
	[0xd32010affc4dab3afc8442132ca50e3e2b62780170f395faacb67ce577b44525, 0xaa7bc42e366f4995960637406ad30445ec0619ad3e898320efdb5a561746129, 0xfe3444472410d50c0538b49086c756139659109b21a64ce2d92b94e591e0ced9, 0x60082c43317965ae18f3bd652447ef5c4dabfd55c54afc0e35054645682ec194, 0x4d4d30e40489b1d453095e65c4290d870c86845ff2eba0f2f246e4b7c6a1c458, 0x7ffd7a904c1602f0a648204e902b0b3d4a46652bb2325e2082b392298c0a2e21, 0x4bbc2aef433c6413eca6c3e6b3b8f0e4d49fd726a9cb659f102fd0fdc0cc58ce, 0x50d7de45d96f38fde92f70f864da20f9596a4371dfb6ae48e5867b09a2ddeefa, 0xe197da8e28cffb508666fd078ec4458faa8e680d77ccf9a27f938a9009aee666, ]
	*/
	
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
		T[i] = Rlo[i];
		T2[i] = Rhi[i];
		
		T[i] += V[i];
		Rlo[i] = V2[i] + T2[i] + (T[i] < V[i]);
		res->lo[i] = LOW(Rlo[i]);
		res->hi[i] = HIGH(Rlo[i]);
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
	
	printf("[");
	for(i = 0; i < N; i++)
		printf("0x%lx%016lx%016lx%016lx, ", HIGH(Rhi[i]), LOW(Rhi[i]), HIGH(Rlo[i]), LOW(Rlo[i]));
	printf("]\n\n");
	
	mns128_montg_int_red(res, Rhi, Rlo);
	
end:
	free_poly(stok);
}


int main(void)
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
	
	int64_t Ahi = 0xffffffffffffffff, Bhi = 0xffffffffffffffff;
	uint64_t  Alo = 0xffffffffffffffff, Blo = 0xffffffffffffffff;
	
	const char RES[] = 
	"0x1b5dc7ca3fcbcc34673dbafa172c2d2ca8fc6b0e76b6a3dd70ff9147383b73";
	
/*	int64_t Ahi = 0x7e4deecd7fb21228, Bhi = -0x11fb30a7f4f269dd;*/
/*	uint64_t  Alo = 0xbe23c803cef47a1, Blo = 0xc97d66ba9ffee37;*/
/*	*/
/*	const char RES[] = */
/*	"0xf720e4b9bc91c7b14365a9592434dee8edbb60fd615d17fd17ab67202e5f1197";*/
	
/*	__int128 R1 = 0;*/
/*	unsigned __int128 R2 = 0;*/
	
	__int128 R1 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;
	unsigned __int128 R2 = (__int128) ((__int128) 0xffffffffffffffff << 64) + 0xffffffffffffffff;
	
	multadd128(&R1, &R2, Ahi, Alo, Bhi, Blo);
	
	printf("%s\n", RES);
	
	printf("0x%lx%016lx%016lx%016lx\n", HIGH(R1), LOW(R1), HIGH(R2), LOW(R2));
	
	printf("\n\n");
	
	poly128 A;
	init_poly128(N, &A);
	const char a[] = "74ff560400d0105e6381e4f7cf22ba4a3d949bbe3b03e7ec1c8aebfb02a4dedf230eef099cd1ae78adf8f142cd70ed93122a5c48c5edcba658615fa2316994dce0c84e9e54c5ae9482acdc0ed6fae84eb7e83d94016d12452ad41369e33a53a676d539439488bdc8b3462c5579a432e8b579e8af9d5b2b0b8f37856fe2de7f30";
	
	convert_string_to_amns128(A, a);
	
	p128_print(A);
	
	free_poly128(A);
	
	return 0;
}

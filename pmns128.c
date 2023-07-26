#include <stdio.h>
#include <time.h>
//#include <gmp.h>
#include <string.h>

#include "pmns128.h"
#include "utilitymp.h"

#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 16")
#define MULTITHREA
#define NBTHREADZ 8

#define AMNS128_MONTG_MULT amns128_montg_mult

static inline void multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const uint64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// Multiplies A and B and adds the result to R using the classic
	// schoolbook algorithm
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux2, aux3;
	
	A1B1 = (__int128) Ahi * Bhi;
	A1B0 = (__int128) Ahi * Blo;
	A0B1 = (__int128) Alo * Bhi;
	A0B0 = (__int128) Alo * Blo;
	
	aux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0);
	
	tmplo = (__int128) A0B0 + ((__int128)(LOW(A0B1) + LOW(A1B0)) << 64);
	*Rhi += (__int128) aux2 + A1B1 + add_overflow(Rlo, tmplo);
}

static inline void mns_multadd128(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const int64_t Ahi, const uint64_t Alo,
	const int64_t Bhi, const uint64_t Blo)
{
	// Multiplies A and B and adds the result to R faster in the case where A
	// and B are both in the PMNS and therefore their coefficients are lower
	// than rho. In that case the middle sum will never overflow, speeding
	// up the calculations.
	unsigned __int128 A0B0, A1B0_A0B1, tmplo;
	__int128 A1B1;
	
	A1B1 = (__int128) Ahi * Bhi;
	A0B0 = (__int128) Alo * Blo;
	A1B0_A0B1 = (__int128) ((__int128) ((__int128) Ahi * Blo) + ((__int128) Alo * Bhi)) + HIGH(A0B0);
	
	tmplo = (__int128) LOW(A0B0) | ((__int128)A1B0_A0B1 << 64);
	*Rhi += (__int128) A1B1 + HIGH(A1B0_A0B1) + add_overflow(Rlo, tmplo);
}

#ifdef LAMBDA

void toepmns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	// Toeplitz version.
	__int128 vect[N], matr[2*N - 1];
	
	const uint8_t CHECKDEX = 0;
	__int128 tmplo, tmphi;
	tmphi = Rhi[CHECKDEX];
	tmplo = Rlo[CHECKDEX];
	
	for(int i = 0; i < N; i++)
	{
		vect[i] = ((__int128)A->hi[i] << 64) | A->lo[i];
		matr[i + N - 1] = ((__int128)B->hi[i] << 64) | B->lo[i];
	}
	for(int i = 1; i < N; i++)
		matr[i - 1] = (__int128) matr[i + N - 1] * LAMBDA;
	
	//toeplitz128_vm3x3_72x72(Rhi, Rlo, vect, matr);
	//toeplitz128_vm3x3_36x36(Rhi, Rlo, vect, matr);
	//toeplitz128_vm3x3_18x18(Rhi, Rlo, vect, matr);
	/*printf("\n\n0x");
	__print128(Rhi[CHECKDEX]); printf(" "); __print128(Rlo[CHECKDEX]);
	Rhi[CHECKDEX] = tmphi; Rlo[CHECKDEX] = tmplo;*/
	
	/*const uint16_t size = 72;
	__int128 A1B0_A0B1;
	register uint16_t i, j;
	
	#pragma omp parallel for num_threads(8) private(A1B0_A0B1)
	for(i = 0; i < size; i++)
	{
		Rhi[i] += (__int128) HIGH(vect[0]) * HIGH(matr[size - 1 + i]);
		for(j = 1; j < size; j++)
			Rhi[i] += (__int128) HIGH(vect[j]) * HIGH(matr[size - 1 - j + i]);
		
		A1B0_A0B1 = ((__int128) HIGH(vect[0]) * LOW(matr[size - 1 + i])) + ((__int128) LOW(vect[0]) * HIGH(matr[size - 1 + i]));
		
		for(j = 1; j < size; j++)
			A1B0_A0B1 += ((__int128) HIGH(vect[j]) * LOW(matr[size - 1 - j + i])) + ((__int128) LOW(vect[j]) * HIGH(matr[size - 1 - j + i]));
		
		Rhi[i] += HIGH(A1B0_A0B1) + add_overflow(Rlo + i, (__int128) ((__int128)A1B0_A0B1 << 64));
		//Rlo[i] += (__int128) ((__int128)A1B0_A0B1 << 64);
		
		for(j = 0; j < size; j++)
			Rhi[i] += add_overflow(Rlo + i , (__int128) LOW(vect[j]) * LOW(matr[size - 1 - j + i]));
	}
		/*printf("\n\n0x");
		__print128(Rhi[CHECKDEX]); printf(" "); __print128(Rlo[CHECKDEX]);
		printf("\n\n"); exit(1);*/
	
}

void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	//register uint16_t i, j;
	unsigned __int128 A1B0_A0B1;
	
	unsigned __int128 A0B0;
	__int128 A1B1, aux1, aux2;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(NBTHREADZ) private(A1B0_A0B1, aux1, aux2)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N - 1; i++)
	{
		A1B0_A0B1 = (__int128) A->hi[i + 1] * B->lo[N - 1] + (__int128) A->lo[i + 1] * B->hi[N - 1];
		Rhi[i] = (__int128) A->hi[i + 1] * B->hi[N - 1];
		Rlo[i] = (__int128) A->lo[i + 1] * B->lo[N - 1];
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 2; j < N - i; j++)
		{
			A1B1 = (__int128) A->hi[i + j] * B->hi[N - j];
			A0B0 = (__int128) A->lo[i + j] * B->lo[N - j];
			A1B0_A0B1 += (__int128) A->hi[i + j] * B->lo[N - j] + (__int128) A->lo[i + j] * B->hi[N - j];
			
			Rhi[i] += (__int128) A1B1 + add_overflow(Rlo + i, A0B0);
		}
		
		Rhi[i] += (__int128) HIGH(A1B0_A0B1) + add_overflow(Rlo + i, (__int128) (LOW(A1B0_A0B1)) << 64);
		
		aux1 = (unsigned __int128) LOW(Rlo[i]) * LAMBDA;
		aux2 = (unsigned __int128) HI(Rlo[i]) * LAMBDA + HIGH(aux1);
		Rlo[i] = ((__int128) aux2 << 64) | LOW(aux1);
		Rhi[i] = (__int128) Rhi[i] * LAMBDA + HIGH(aux2);
		
		
		A1B0_A0B1 = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			A1B1 = (__int128) A->hi[j] * B->hi[i - j];
			A0B0 = (__int128) A->lo[j] * B->lo[i - j];
			A1B0_A0B1 += (__int128) ((__int128) ((__int128) A->hi[j] * B->lo[i - j]) + ((__int128) A->lo[j] * B->hi[i - j]));
	
			Rhi[i] += (__int128) A1B1 + add_overflow(Rlo + i, A0B0);
		}
		
		Rhi[i] += (__int128) HIGH(A1B0_A0B1) + add_overflow(Rlo + i, (__int128) (LOW(A1B0_A0B1)) << 64);
		
	}
	A1B0_A0B1 = (__int128) ((__int128) ((__int128) A->hi[0] * B->lo[N - 1]) + ((__int128) A->lo[0] * B->hi[N - 1]));
	Rhi[N - 1] = (__int128) A->hi[0] * B->hi[N - 1];
	Rlo[N - 1] = (__int128) A->lo[0] * B->lo[N - 1];
	_PRAGMAGCCUNROLLLOOP_
	for(int j = 1; j < N; j++)
	{
		A1B1 = (__int128) A->hi[j] * B->hi[N - 1 - j];
		A0B0 = (__int128) A->lo[j] * B->lo[N - 1 - j];
		A1B0_A0B1 += (__int128) ((__int128) ((__int128) A->hi[j] * B->lo[N - 1 - j]) + ((__int128) A->lo[j] * B->hi[N - 1 - j]));
	
		Rhi[N - 1] += (__int128) A1B1 + add_overflow(Rlo + N - 1, A0B0);
	}
	
	Rhi[N - 1] += (__int128) HIGH(A1B0_A0B1) + add_overflow(Rlo + N - 1, (__int128) (LOW(A1B0_A0B1)) << 64);
	
}

#endif

#ifdef LENEXTPOLY

void mns128_mod_mult_ext_red(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) an irreducible polynomial used for reduction. Result in R
	__int128 Thi;
	unsigned __int128 aux1, aux2, Tlo;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(8) private(Thi, Tlo, aux1, aux2)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		Thi = 0;
		Tlo = 0;
		for(int j = 1; j < N - i; j++)
			mns_multadd128(&Thi, &Tlo, A->hi[i + j], A->lo[i + j],
				B->hi[N - j], B->lo[N - j]);
		for(int k = 0; (k < LENEXTPOLY) && (i + k < N); k++)
		{
			aux1 = (unsigned __int128) LOW(Tlo) * EXTPOLY[k] + LOW(Rlo[i + k]);
			aux2 = (unsigned __int128) HI(Tlo) * EXTPOLY[k] + HIGH(aux1) + HI(Rlo[i + k]);
			Rlo[i + k] = ((__int128) aux2 << 64) | LOW(aux1);
			Rhi[i + k] += (__int128) Thi * EXTPOLY[k] + HIGH(aux2);
		}
		
		for(int j = 0; j < i + 1; j++)
			mns_multadd128(Rhi + i, Rlo + i, A->hi[j], A->lo[j],
				B->hi[i - j], B->lo[i - j]);
	}
}

#endif

#ifdef M_or_B_is_M

void m_or_b_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	//register uint16_t i, j;
	
	//#pragma omp parallel for num_threads(8)
	/*for(int i = 0; i < N; i++)
	{
		Rhi[i] += Rlo[i] != 0;
		Rlo[i] = 0;
		for(int j = 1; j < N - i; j++)
			multadd128(Rhi + i, Rlo + i, HIGH(A[i + j]), LOW(A[i + j]),
				MLambdahi[N - j], MLambdalo[N - j]);
		
		for(int j = 0; j < i + 1; j++)
			multadd128(Rhi + i, Rlo + i, HIGH(A[j]), LOW(A[j]),
				Mhi[i - j], Mlo[i - j]);
	}*/
	
	//register uint16_t i, j;
	unsigned __int128 A0B1;
	__int128 A1B0;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(NBTHREADZ) private(A0B1, A1B0)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			Rhi[i] += (__int128) HIGH(A[j]) * HIGH(matrM[N - 1 - j + i]) + add_overflow(Rlo + i, ((__int128) LOW(A[j]) * LOW(matrM[N - 1 - j + i])));
		
		A0B1 = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			A0B1 += ((__int128) LOW(A[j]) * HIGH(matrM[N - 1 - j + i]));
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
		{
			A1B0 = (__int128) HIGH(A[j]) * LOW(matrM[N - 1 - j + i]);
			Rhi[i] += HIGH(A1B0) + add_overflow(Rlo + i, ((__int128)(LOW(A1B0)) << 64));
		}
		
		
		Rhi[i] += HIGH(A0B1) + add_overflow(Rlo + i, ((__int128)LOW(A0B1) << 64));
	}
}

void toepm_or_b_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	/*const uint8_t CHECKDEX = 24;
	__int128 tmplo, tmphi;
	tmphi = Rhi[CHECKDEX];
	tmplo = Rlo[CHECKDEX];*/
	//Mtoeplitz128_vm3x3_18x18(Rhi, Rlo, A, matrM);
	
	//Mtoeplitz128_vm3x3_72x72(Rhi, Rlo, A, matrM);
	//Mtoeplitz128_vm3x3_36x36(Rhi, Rlo, A, matrM);
	//Mtoeplitz128_vm(Rhi, Rlo, A, matrM);
	
	/*printf("\n\n0x");
	__print128(Rhi[CHECKDEX]); printf(" "); __print128(Rlo[CHECKDEX]);
	Rhi[CHECKDEX] = tmphi; Rlo[CHECKDEX] = tmplo;
	
	normm_or_b_mns128_mod_mult_ext_red(Rhi, Rlo, A);
	printf("\n\n0x");
	__print128(Rhi[CHECKDEX]); printf(" "); __print128(Rlo[CHECKDEX]);
	printf("\n\n"); exit(1);*/
}

void m1_or_b1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R. 
	// In pmns128 we only care about the lower 128 bits for this operation.
	
	//register uint16_t i, j;
	
	//#pragma omp parallel for num_threads(8)
	/*for(int i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		for(int j = 1; j < N - i; j++)
			Rlo[i] += (__int128) A[i + j] * M1Lambda[N - j]; //(((__int128) LOW(B1hi[j][i]) << 64) | B1lo[j][i]));
			//m1_multadd128(Rlo + i, A[i+j], M1Lambda[N - j]);
		
		for(int j = 0; j < i + 1; j++)
			Rlo[i] += (__int128) A[j] * M1[i - j];
			//m1_multadd128(Rlo + i, A[j], M1[i - j]);
	}*/
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		Rlo[i] = 0;
		
		for(int j = 0; j < N; j++)
			Rlo[i] += (__int128) A[j] * matrM1[N - 1 - j + i];
	}
}

void toepm1_or_b1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R. 
	// In pmns128 we only care about the lower 128 bits for this operation.
	// Toeplitz version.

	//m1toeplitz128_vm3x3_18x18(Rlo, A, matrM1);
	//m1toeplitz128_vm3x3_72x72(Rlo, A, matrM1);
	//m1toeplitz128_vm3x3_36x36(Rlo, A, matrM1);
}


#endif

#ifdef M_or_B_is_B

void m_or_b_mns128_mod_mult_ext_red(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, unsigned __int128* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	
	unsigned __int128 A0B1;
	__int128 A1B0, aux;
	
	//_PRAGMAGCCUNROLLLOOP_
	#ifdef MULTITHREAD
	#pragma omp parallel for private(A0B1, A1B0, aux) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		aux = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			Rhi[i] += (__int128) HIGH(A[j]) * Bhi[j][i] + add_overflow(Rlo + i, ((__int128) LOW(A[j]) * Blo[j][i]));
		
		A0B1 = (__int128) LOW(A[0]) * Bhi[0][i];
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N; j++)
			A0B1 += (__int128) LOW(A[j]) * Bhi[j][i];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
		{
			A1B0 = (__int128) HIGH(A[j]) * Blo[j][i];
/*			aux += add_overflow(&A0B1, A1B0) * (1 - 2 * (A1B0 < 0));*/
			Rhi[i] += HIGH(A1B0) + add_overflow(Rlo + i, ((__int128)(LOW(A1B0)) << 64));
		}
		
/*		A0B1 = 0;*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = 0; j < N; j++)*/
/*			A0B1 += ((__int128) LOW(A[j]) * Bhi[j][i]);*/
		
		Rhi[i] += (__int128) HIGH(A0B1) + add_overflow(Rlo + i, ((__int128)LOW(A0B1) << 64));
		
/*		Rhi[i] += Rlo[i] != 0;*/
/*		Rlo[i] = 0;*/
/*		*/
/*		_PRAGMAGCCUNROLLLOOP_*/
/*		for(j = 0; j < N; j++)*/
/*			multadd128(Rhi + i, Rlo + i, HIGH(A[j]), LOW(A[j]),*/
/*				Bhi[j][i], Blo[j][i]);*/
	}
}

void m1_or_b1_mns128_mod_mult_ext_red(unsigned __int128* restrict Rlo,
	unsigned __int128* restrict A)
{
	// Vector-Matrix multiplication between A and B1, result in R.
	// In 128 bit version we only care about the lower 128 bits.
	
	
	#ifdef MULTITHREAD
	#pragma omp parallel for num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		Rlo[i] = ((unsigned __int128) A[0] * (((unsigned __int128) B1hi[0][i] << 64) | B1lo[0][i]));
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N; j++)
			Rlo[i] += ((unsigned __int128) A[j] * (((unsigned __int128) B1hi[j][i] << 64) | B1lo[j][i]));
	}
}

#endif

void mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Function that reduces the internal coefficient contained in R to be lower
	// than the chosen Rho using Montgomery's internal reduction algorithm.
	unsigned __int128 V[N];
	register uint16_t i;
	
	m1_or_b1_mns128_mod_mult_ext_red(V, Rlo);
	
	m_or_b_mns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}

/*void UNROLLED_mns128_montg_int_red(poly128 res, __int128* restrict Rhi,
	unsigned __int128* restrict Rlo)
{
	// Unrolled version.
	unsigned __int128 V[N];
	register uint16_t i;
	
	//UNROLLED_m1_or_b1_mns128_mod_mult_ext_red(V, Rlo);
	
	//UNROLLED_m_or_b_mns128_mod_mult_ext_red(Rhi, Rlo, V);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Rhi[i]);
		res->hi[i] = HIGH(Rhi[i]);
	}
}*/

void Xamns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	mns128_montg_int_red(res, Rhi, Rlo);
}

#ifdef M_or_B_is_M

void toeplitz_vm_24x24(int64_t* rophi, __int128* restrict rop, const uint64_t* restrict vect, const __int128* restrict matr)
{
	for(int i = 0; i < 24; i++)
	{
		rop[i] = (__int128) vect[0] * matr[24 - 1 + i];
		rophi[i] = __builtin_add_overflow(rop[i], (__int128) vect[1] * matr[24 - 2 + i], rop + i);
		for(int j = 2; j < 24; j++)
			rophi[i] += __builtin_add_overflow(rop[i], (__int128) vect[j] * matr[24 - 1 - j + i], rop + i);
	}
}

void utoeplitz_vm_24x24(int64_t* rophi, __int128* restrict rop, const __int128* restrict vect, const __int128* restrict matr)
{
	for(int i = 0; i < 24; i++)
	{
		rop[i] = (__int128) vect[0] * matr[24 - 1 + i];
		rophi[i] = __builtin_add_overflow(rop[i], (__int128) vect[1] * matr[24 - 2 + i], rop + i);
		for(int j = 2; j < 24; j++)
			rophi[i] += __builtin_add_overflow(rop[i], (__int128) vect[j] * matr[24 - 1 - j + i], rop + i);
	}
}

void toeplitz128_vm_24x24(uint64_t* restrict rophi, __int128* restrict roplo, const __int128* restrict vect, const __int128* restrict matr)
{
	/*__int128 t0hi[12], t1hi[12], t2hi[12];
	unsigned __int128 t0lo[12], t1lo[12], t2lo[12];
	__int128 v0p1[12], m0m1[23], m0m2[23];
	register uint16_t i;
	
	for(i = 0; i < 12; i++)
	{
		v0p1[i] = vect[i] + vect[i + 12];
	}
	
	for(i = 0; i < 23; i++)
	{
		m0m1[i] = matr[12 + i] - matr[24 + i];
		m0m2[i] = matr[12 + i] - matr[i];
	}
	
	toeplitz128_vm_12x12(t0hi, t0lo, v0p1, matr + 12);
	toeplitz128_vm_12x12(t1hi, t1lo, vect, m0m1);
	toeplitz128_vm_12x12(t2hi, t2lo, vect + 12, m0m2);
	
	for(i = 0; i < 12; i++)
	{
		rophi[i] = t0hi[i] - t2hi[i] - __builtin_sub_overflow(t0lo[i], t2lo[i], roplo + i);
		rophi[i + 12] = t0hi[i] - t1hi[i] - __builtin_sub_overflow(t0lo[i], t1lo[i], roplo + i + 12);
	}*/
	for(int i = 0; i < 24; i++)
	{
		roplo[i] = 0;
		rophi[i] = 0;
		for(int j = 0; j < 24; j++)
		{
			rophi[i] += __builtin_add_overflow((__int128)roplo[i], (__int128) vect[j] * matr[24 - 1 - j + i], roplo + i);
		}
	}
}

void toeplitz128_vm3x3_72x72(uint64_t* restrict rophi, __int128* restrict roplo, const uint64_t* restrict vect, const __int128* restrict matr)
{
	uint64_t t0hi[24], t1hi[24], t2hi[24], t3hi[24], t4hi[24], t5hi[24];
	__int128 t0lo[24], t1lo[24], t2lo[24], t3lo[24], t4lo[24], t5lo[24];
	
	__int128 m03, v0[24], v1[24], v2[24], v1m2[24], v0m2[24], v0m1[24], m034[47], m013[47], m012[47], m0[47], m1[47], m3[47];
	register uint16_t i;
	
	for(i = 0; i < 24; i++)
	{
		v0[i] = vect[i];
		v1[i] = vect[24 + i];
		v2[i] = vect[48 + i];
		v1m2[i] = (__int128) vect[24 + i] - vect[48 + i];
		v0m1[i] = (__int128) vect[i] - vect[24 + i];
		v0m2[i] = (__int128) vect[i] - vect[48 + i];
	}
	
	for(i = 0; i < 47; i++)
	{
		m0[i] = matr[i + 48];
		m1[i] = matr[i + 72];
		m3[i] = matr[i + 24];
		m03 = (__int128) m0[i] + m3[i];
		m034[i] = (__int128) m03 + matr[i];
		m013[i] = (__int128) m03 + m1[i];
		m012[i] = (__int128) m0[i] + m1[i] + matr[i + 96];
	}
	
	toeplitz128_vm_24x24(t0hi, t0lo, v2, m034);
	toeplitz128_vm_24x24(t1hi, t1lo, v1, m013);
	toeplitz128_vm_24x24(t2hi, t2lo, v0, m012);
	toeplitz128_vm_24x24(t3hi, t3lo, v1m2, m3);
	toeplitz128_vm_24x24(t4hi, t4lo, v0m2, m0);
	toeplitz128_vm_24x24(t5hi, t5lo, v0m1, m1);
	
	for(i = 0; i < 24; i++)
	{
		rophi[i] = t0hi[i] + t3hi[i] + t4hi[i] + __builtin_add_overflow(t0lo[i], t3lo[i], roplo + i) + __builtin_add_overflow(roplo[i], t4lo[i], roplo + i);
		rophi[i + 24] = t1hi[i] - t3hi[i] + t5hi[i] - __builtin_sub_overflow(t1lo[i], t3lo[i], roplo + i + 24) + add_overflow(roplo + i + 24, t5lo[i]);
		rophi[i + 48] = t2hi[i] - t5hi[i] - t4hi[i] - __builtin_sub_overflow(t2lo[i], t4lo[i], roplo + i + 48) - __builtin_sub_overflow(roplo[i + 48], t5lo[i], roplo + i + 48);
	}
}

void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the Montgomery CIOS method in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will be in the PMNS also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	int64_t tmplo;
	uint64_t Shi[N];
	unsigned __int128 tmp, aux, Tlo[N] = {0}, Thi[N] = {0};
	__int128 Slo[N];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N - 1; i++)
	{
		Slo[i] = (unsigned __int128) A->lo[i + 1] * B->lo[N - 1];
		Shi[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 2; j < N - i; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[i + j] * B->lo[N - j]);
		}
		
		aux = (unsigned __int128) LOW(Slo[i]) * LAMBDA;
		tmp = (unsigned __int128) HI(Slo[i]) * LAMBDA + HI(aux);
		Slo[i] = ((unsigned __int128) tmp << 64) | LOW(aux);
		Shi[i] = (unsigned __int128) Shi[i] * LAMBDA + HI(tmp);
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[j] * B->lo[i - j]);
		}
	}
	Slo[N - 1] = (unsigned __int128) A->lo[0] * B->lo[N - 1];
	Shi[N - 1] = add_overflow(Slo + N - 1, (unsigned __int128) A->lo[1] * B->lo[N - 2]);
	_PRAGMAGCCUNROLLLOOP_
	for(int j = 2; j < N; j++)
	{
		Shi[N - 1] += add_overflow(Slo + N - 1, (unsigned __int128) A->lo[j] * B->lo[N - 1 - j]);
	}
	
	printf("(0x%lu << 128) | (0x", Shi[0]); __print128(Slo[0]); printf(")\n");
	
	/*{
		uint64_t *vect = A->lo;
		__int128 matr[143];
		
		for(int i = 0; i < 71; i++)
		{
			matr[i + 71] = B->lo[i];
			matr[i] = B->lo[1 + i] * LAMBDA;
		}
		matr[142] = B->lo[71];
		
		int64_t t0hi[24], t1hi[24], t2hi[24], t3hi[24], t4hi[24], t5hi[24];
		__int128 t0[24], t1[24], t2[24], t3[24], t4[24], t5[24];
		
		__int128 m03, v1m2[24], v0m2[24], v0m1[24], m034[47], m013[47], m012[47], m0[47], m1[47], m3[47];
		
		for(int i = 0; i < 24; i++)
		{
			v1m2[i] = (__int128) vect[24 + i] - vect[48 + i];
			v0m1[i] = (__int128) vect[i] - vect[24 + i];
			v0m2[i] = (__int128) vect[i] - vect[48 + i];
		}
		
		for(int i = 0; i < 47; i++)
		{
			m0[i] = matr[i + 48];
			m1[i] = matr[i + 72];
			m3[i] = matr[i + 24];
			m03 = (__int128) m0[i] + m3[i];
			m034[i] = (__int128) m03 + matr[i];
			m013[i] = (__int128) m03 + m1[i];
			m012[i] = (__int128) m0[i] + m1[i] + matr[i + 96];
		}
		
		toeplitz_vm_24x24(t0hi, t0, vect + 56, m034);
		toeplitz_vm_24x24(t1hi, t1, vect + 28, m013);
		toeplitz_vm_24x24(t2hi, t2, vect, m012);
		utoeplitz_vm_24x24(t3hi, t3, v1m2, m3);
		utoeplitz_vm_24x24(t4hi, t4, v0m2, m0);
		utoeplitz_vm_24x24(t5hi, t5, v0m1, m1);
		
		for(int i = 0; i < 24; i++)
		{
			Shi[i] = t0hi[i] + t3hi[i] + t4hi[i] + __builtin_add_overflow(t0[i], t3[i], Slo + i) + __builtin_add_overflow(Slo[i], t4[i], Slo + i);
			Shi[i + 24] =  t1hi[i] - t3hi[i] + t5hi[i] + __builtin_add_overflow(t1[i], t5[i], Slo + i + 24) - __builtin_sub_overflow(Slo[i + 24], t3[i], Slo + i + 24);
			Shi[i + 48] =  t2hi[i] - t4hi[i] - t5hi[i] - __builtin_sub_overflow(t2[i], t4[i], Slo + i + 48) - __builtin_sub_overflow(Slo[i + 48], t5[i], Slo + i + 48);
		}
		for(i = 0; i < 24; i++)
	{
		rophi[i] = t0hi[i] + t3hi[i] + t4hi[i] + __builtin_add_overflow(t0lo[i], t3lo[i], roplo + i) + add_overflow(roplo + i, t4lo[i]);
		rophi[i + 24] = t1hi[i] - t3hi[i] + t5hi[i] - __builtin_sub_overflow(t1lo[i], t3lo[i], roplo + i + 24) + add_overflow(roplo + i + 24, t5lo[i]);
		rophi[i + 48] = t2hi[i] - t5hi[i] - t4hi[i] - __builtin_sub_overflow(t2lo[i], t4lo[i], roplo + i + 48) - sub_overflow(roplo + i + 48, t5lo[i]);
	}
	}*/
	
	__int128 matr[143];
	
	for(int i = 0; i < 71; i++)
	{
		matr[i + 71] = B->lo[i];
		matr[i] = B->lo[1 + i] * LAMBDA;
	}
	matr[142] = B->lo[71];
	toeplitz128_vm3x3_72x72(Shi, Slo, A->lo, matr);
	
	
	printf("(0x%lu << 128) | (0x", Shi[0]);
	__print128(Slo[0]);
	printf(")\n");
	exit(0);
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * matrM1[N - 1 + i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * matrM1[N - 1 - j + i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) LOW(matrM[N - 1 - i + j]) * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) HIGH(matrM[N - 1 - i + j]) * tmplo + HIGH(tmp);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HI(Slo[i]) + ((__int128)Shi[i] << 64) + HIGH(Tlo[i]);
		Slo[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
		{
			Slo[i] += (__int128) A->lo[i + j] * B->hi[N - j] +
				(__int128) A->hi[i + j] * B->lo[N - j];
		}
		
		Slo[i] *= LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			Slo[i] += (__int128) A->hi[j] * B->lo[i - j] +
				(__int128) A->lo[j] * B->hi[i - j];
		}
		Slo[i] += Thi[i] + 1;
		Thi[i] = 0;
		Tlo[i] = 0;
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * matrM1[N - 1 + i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * matrM1[N - 1 - j + i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) LOW(matrM[N - 1 - i + j]) * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) HIGH(matrM[N - 1 - i + j]) * tmplo + HIGH(tmp);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HIGH(Slo[i]) + HIGH(Tlo[i]) + 1;
		tmp = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
		{
			tmp += (__int128) A->hi[i + j] * B->hi[N - j];
		}
		
		tmp *= LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			tmp += (__int128) A->hi[j] * B->hi[i - j];
		}
		Thi[i] += tmp;
	}
	
	for(int i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Thi[i]);
		res->hi[i] = HIGH(Thi[i]);
	}
}

#endif

#ifdef M_or_B_is_B

void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Function that multiplies A by B using the Montgomery CIOS method in a
	// PMNS. Puts the result in res. A and B have to be in the system and res
	// will be in the PMNS also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	int64_t tmplo;
	uint64_t Shi[N];
	unsigned __int128 tmp, aux, Slo[N], Tlo[N] = {0}, Thi[N] = {0};
	
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N - 1; i++)
	{
		Slo[i] = (unsigned __int128) A->lo[i + 1] * B->lo[N - 1];
		Shi[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 2; j < N - i; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[i + j] * B->lo[N - j]);
		}
		
		aux = (unsigned __int128) LOW(Slo[i]) * LAMBDA;
		tmp = (unsigned __int128) HI(Slo[i]) * LAMBDA + HI(aux);
		Slo[i] = ((unsigned __int128) tmp << 64) | LOW(aux);
		Shi[i] = (unsigned __int128) Shi[i] * LAMBDA + HI(tmp);
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			Shi[i] += add_overflow(Slo + i, (unsigned __int128) A->lo[j] * B->lo[i - j]);
		}
	}
	Slo[N - 1] = (unsigned __int128) A->lo[0] * B->lo[N - 1];
	Shi[N - 1] = add_overflow(Slo + N - 1, (unsigned __int128) A->lo[1] * B->lo[N - 2]);
	_PRAGMAGCCUNROLLLOOP_
	for(int j = 2; j < N; j++)
	{
		Shi[N - 1] += add_overflow(Slo + N - 1, (unsigned __int128) A->lo[j] * B->lo[N - 1 - j]);
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * B1lo[0][i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * B1lo[j][i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[i][j] * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) Bhi[i][j] * tmplo + HIGH(tmp);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HI(Slo[i]) + ((__int128)Shi[i] << 64) + HIGH(Tlo[i]);
		Slo[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
		{
			Slo[i] += (__int128) A->lo[i + j] * B->hi[N - j] +
				(__int128) A->hi[i + j] * B->lo[N - j];
		}
		
		Slo[i] *= LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			Slo[i] += (__int128) A->hi[j] * B->lo[i - j] +
				(__int128) A->lo[j] * B->hi[i - j];
		}
		Slo[i] += Thi[i] + 1;
		Thi[i] = 0;
		Tlo[i] = 0;
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		tmplo = Slo[0] * B1lo[0][i];
		for(int j = 1; j < N; j++)
			tmplo += Slo[j] * B1lo[j][i];
		
		for(int j = 0; j < N; j++)
		{
			tmp = (unsigned __int128) Blo[i][j] * tmplo;
			Tlo[j] += LOW(tmp);
			Thi[j] += (__int128) Bhi[i][j] * tmplo + HIGH(tmp);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N; i++)
	{
		Thi[i] += HIGH(Slo[i]) + HIGH(Tlo[i]) + 1;
		tmp = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
		{
			tmp += (__int128) A->hi[i + j] * B->hi[N - j];
		}
		
		tmp *= LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
		{
			tmp += (__int128) A->hi[j] * B->hi[i - j];
		}
		Thi[i] += tmp;
	}
	
	for(int i = 0; i < N; i++)
	{
		res->lo[i] = LOW(Thi[i]);
		res->hi[i] = HIGH(Thi[i]);
	}
}

#endif

/*void UNROLLED_amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B)
{
	// Unrolled version.
	
	__int128 Rhi[N] = {0};
	unsigned __int128 Rlo[N] = {0};
	
	//UNROLLED_mns128_mod_mult_ext_red(Rhi, Rlo, A, B);
	
	//UNROLLED_mns128_montg_int_red(res, Rhi, Rlo);
}*/

void amns128_rtl_sqandmult(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Right to Left version. Slower than ltr.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly128 tmp;
	init_poly128(N, &tmp);
	
	p128_copy(res, &__theta__);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		for(j = 0; j < 64; j++)
			{
				if(aux & (1ULL << j))
					AMNS128_MONTG_MULT(tmp, res, base);
				else
					AMNS128_MONTG_MULT(tmp, res, &__theta__);
				AMNS128_MONTG_MULT(res, tmp, tmp);
			}
	}
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			AMNS128_MONTG_MULT(tmp, res, base);
		else
			AMNS128_MONTG_MULT(tmp, res, &__theta__);
		AMNS128_MONTG_MULT(res, tmp, tmp);
		aux >>= 1;
	}
	
	free_poly128(tmp);
}

void amns128_ltr_sqandmult(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Left to Right version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly128 tmp;
	init_poly128(N, &tmp);
	
	p128_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		AMNS128_MONTG_MULT(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			AMNS128_MONTG_MULT(res, tmp, base);
		else
			AMNS128_MONTG_MULT(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			AMNS128_MONTG_MULT(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				AMNS128_MONTG_MULT(res, tmp, base);
			else
				AMNS128_MONTG_MULT(res, tmp, &__theta__);
		}
	}
	
	free_poly128(tmp);
}

void amns128_montg_ladder(restrict poly128 res, const restrict poly128 base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux, b;
	
	poly128 tmp, R[2];
	init_poly128(N, &tmp);
	R[0] = res;
	R[1] = tmp;
	
	p128_copy(res, &__theta__);
	p128_copy(tmp, base);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		AMNS128_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
		AMNS128_MONTG_MULT(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			AMNS128_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
			AMNS128_MONTG_MULT(R[b], R[b], R[b]);
		}
	}
	
	free_poly128(tmp);
}

void amns128_sqandmult(restrict poly128 res, const char* restrict base,
	const char* restrict exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Uses montgomery ladder.
	
	poly128 pbase;
	mpnum mpexponent;
	init_poly128(N, &pbase);
	init_mpnum(1, &mpexponent);
	
	convert_string_to_amns128(pbase, base);
	convert_string_to_multipre(&mpexponent, exponent);
	
	amns128_montg_ladder(res, pbase, mpexponent);
	
	free_poly128(pbase);
	free_mpnum(mpexponent);
}

void convert_string_to_amns128(restrict poly128 res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (we multiply it by PHI in the
	// process).
	
	uint8_t counter;
	register uint16_t i, j;
	const unsigned __int128 rho = ((__int128)1) << RHO;
	unsigned __int128 limb, Rlo[N] = {0};
	__int128 Rhi[N] = {0}, tmp[N] = {0};
	mpnum stok, aux1, aux2;
	init_mpnums(2 * N, &stok, &aux1, &aux2, NULL);
	
	convert_string_to_multipre(&aux1, string);
	
	mp_mod(&aux2, aux1, &__P__);
	
	for(i = 0; i < aux2->deg; i++)
		stok->t[i] = aux2->t[i];
	
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
	
	free_mpnums(stok, aux1, aux2, NULL);
}

void convert_amns128_to_multipre(restrict mpnum* res, const restrict poly128 P)
{
	// Function that converts out of the AMNS system and into a multiprecision
	// number.
	
	register uint16_t i;
	mpnum aux, ag, tmp, tmp2;
	poly128 a;
	
	init_mpnums(2 * N, &ag, &tmp, &tmp2, NULL);
	init_mpnum(2, &aux);
	init_poly128(N, &a);
	if((*res)->deg < 2 * N)
	{
		free_mpnum(*res);
		init_mpnum(2 * N, res);
	}
	(*res)->deg = 2 * N;
	for(i = 1; i < 2 * N; i++)
		(*res)->t[i] = 0;
	
	poly128 one;
	init_poly128(N, &one);
	for(i = 1; i < N; i++)
	{
		one->lo[i] = 0;
		one->hi[i] = 0;
	}
	one->hi[0] = 0;
	one->lo[0] = 1;
	
	AMNS128_MONTG_MULT(a, P, one);
	
	free_poly128(one);
	
	(*res)->sign = 1 - 2 * (a->hi[0] < 0);
	(*res)->t[0] = a->lo[0] * (*res)->sign;
	(*res)->t[1] = a->hi[0] * (*res)->sign - ((*res)->sign == -1);
	for(i = 1; i < N; i++)
	{
		aux->sign = 1 - 2 * (a->hi[i] < 0);
		aux->t[0] = a->lo[i] * aux->sign;
		aux->t[1] = a->hi[i] * aux->sign - (aux->sign == -1);
		
		mp_mult(&ag, aux, &Gi[i - 1]);
		
		mp_copy(&tmp, *res);
		mp_add(res, tmp, ag);
	}
	mp_copy(&tmp, *res);
	mp_mod(res, tmp, &__P__);
	
	free_mpnums(aux, ag, tmp, tmp2, NULL);
	free_poly128(a);
}

static inline int64_t randomint64(void)
{
	// Function to generate a random 64 bit number.
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrhohi(int64_t param)
{
	// Utility function to get a high part for a random poly128.
	if (RHO <= 64)
		return 0;
	else
		return param & ((1ULL << (RHO - 64)) - 1);
}

void randpoly128(poly128 P)
{
	// Generates a random poly128 with appropriate, lower than rho high part.
	for(register uint16_t i = 0; i < P->deg; i++)
	{
		P->lo[i] = randomint64();
		if(RHO > 64) P->hi[i] = __modrhohi(randomint64()) * (1 + (rand() & 1) * -2);
	}
}

void __multchecks__(char* nbmults)
{
	// Debug tool to quickly check if the pmns gives the correct results.
	poly128 a, b, c;
	init_poly128s(N, &a, &b, &c, NULL);
	int64_t seed;
	
	register uint64_t cap = 100;
	
	if(nbmults[0] != '\0')
		cap = atoll(nbmults);
	
	//srand((unsigned) (time(&seed)));
	
	for(register uint64_t i = 0; i < cap; i++)
	{
		randpoly128(a);
		randpoly128(b);
		p128_print(a);
		p128_print(b);
		AMNS128_MONTG_MULT(c, a, b);
		p128_print(c);
	}
	
	free_poly128s(a, b, c, NULL);
}

void __sqandmultdemo(void)
{
	// TODO: delete it, used to check point to point process.
	const char a[] = "0xffffffffffffffffffffffffffffffffffffffffffffffff",
		b[] = "0xffffffffffffffffffffffffffffffffffffffffffffffff",
		c[] = "0xfca02b1ec94f1d1a5afdbacd6c0b04b9ded8463b62fade4d35e64f3de3a20e8b172382d318083548be85676b8b7dd347c32841c3efa88c928e7cc945f5e2253708053c8f75ff222717ea374b7f0b6ceafd6ac0b119d23300ac9b6a23cd0b332b438f7be9f6c4ee4cdbbe93444343a8011583e10a8c4cc0f07b67d1467b9a73e";
	
	poly128 C;
	mpnum aux;
	init_poly128(N, &C);
	init_mpnum(1, &aux);
	
	amns128_sqandmult(C, a, b);
	
	convert_amns128_to_multipre(&aux, C);
	
	printf("\n%s\n\n", c);
	mp_print(aux);
	printf("\n");
	
	free_poly128(C);
	free_mpnum(aux);
}

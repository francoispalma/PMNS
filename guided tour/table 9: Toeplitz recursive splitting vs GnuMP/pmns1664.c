#include <stdint.h>

#include "pmns1664.h"
#include "params1664.h"

const mpnum __P1664__ = &__P__;
const mpnum Gamma1664 = Gi;
const int8_t N1664 = N;
const int8_t RHO1664 = RHO;

static inline void toeplitz_vm_8x8(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 8
	for(int i = 0; i < 8; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 8
		for(int j = 0; j < 8; j++)
			rop[i] += (__int128) vect[j] * matr[8 - 1 - j + i];
	}
}

static inline void toeplitz_vm_16x16(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[8], t1[8], t2[8];
	int64_t v0p1[8], m0m1[15], m0m2[15];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
		v0p1[i] = vect[i] + vect[i + 8];
	#pragma GCC unroll 96
	for(int i = 0; i < 15; i++)
	{
		m0m1[i] = matr[i + 8] - matr[i + 16];
		m0m2[i] = matr[i + 8] - matr[i];
	}
	
	toeplitz_vm_8x8(t0, v0p1, matr + 8);
	toeplitz_vm_8x8(t1, vect, m0m1);
	toeplitz_vm_8x8(t2, vect + 8, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 8] = t0[i] - t1[i];
	}
}

static inline void toeplitz_vm_32x32(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[16], t1[16], t2[16];
	int64_t v0p1[16], m0m1[31], m0m2[31];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
		v0p1[i] = vect[i] + vect[i + 16];
	#pragma GCC unroll 96
	for(int i = 0; i < 31; i++)
	{
		m0m1[i] = matr[i + 16] - matr[i + 32];
		m0m2[i] = matr[i + 16] - matr[i];
	}
	
	toeplitz_vm_16x16(t0, v0p1, matr + 16);
	toeplitz_vm_16x16(t1, vect, m0m1);
	toeplitz_vm_16x16(t2, vect + 16, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 16] = t0[i] - t1[i];
	}
}

static inline void Movtoeplitz_vm_4x4(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 4
	for(int i = 0; i < 4; i++)
	{
		#pragma GCC unroll 4
		for(int j = 0; j < 4; j++)
			rop[i] += (__int128) vect[j] * matr[4 - 1 - j + i];
	}
}

static inline void Movtoeplitz_vm_8x8(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[4] = {0};
	int64_t m0m1[7], m0m2[7];
	__int128 v0p1[4];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 4; i++)
	{
		v0p1[i] = (__int128) vect[i] + vect[i + 4];
	}
	#pragma GCC unroll 96
	for(int i = 0; i < 7; i++)
	{
		m0m1[i] = -matr[i + 4] + matr[i + 8];
		m0m2[i] = -matr[i + 4] + matr[i];
	}
	
	Movtoeplitz_vm_4x4(t0, v0p1, matr + 4);
	Movtoeplitz_vm_4x4(rop + 4, vect, m0m1);
	Movtoeplitz_vm_4x4(rop, vect + 4, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 4; i++)
	{
		rop[i] += t0[i];
		rop[i + 4] += t0[i];
	}
}

static inline void Movtoeplitz_vm_16x16(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[8] = {0};
	int64_t m0m1[15], m0m2[15];
	__int128 v0p1[8];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
		v0p1[i] = (__int128) vect[i] + vect[i + 8];
	#pragma GCC unroll 96
	for(int i = 0; i < 15; i++)
	{
		m0m1[i] = -matr[i + 8] + matr[i + 16];
		m0m2[i] = -matr[i + 8] + matr[i];
	}
	
	Movtoeplitz_vm_8x8(t0, v0p1, matr + 8);
	Movtoeplitz_vm_8x8(rop + 8, vect, m0m1);
	Movtoeplitz_vm_8x8(rop, vect + 8, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
	{
		rop[i] += t0[i];
		rop[i + 8] += t0[i];
	}
}

static inline void Mtoeplitz_vm_16x16(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 16
	for(int i = 0; i < 16; i++)
	{
		#pragma GCC unroll 16
		for(int j = 0; j < 16; j++)
			rop[i] += (__int128) vect[j] * matr[16 - 1 - j + i];
	}
}

static inline void Mtoeplitz_vm_32x32(__int128* restrict rop, const int64_t* restrict vect)
{
	__int128 t0[16] = {0};
	int64_t m0m1[31], m0m2[31];
	__int128 v0p1[16];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
		v0p1[i] = (__int128) vect[i] + vect[i + 16];
	#pragma GCC unroll 96
	for(int i = 0; i < 31; i++)
	{
		m0m1[i] = -matrM[i + 16] + matrM[i + 32];
		m0m2[i] = -matrM[i + 16] + matrM[i];
	}
	
	Movtoeplitz_vm_16x16(t0, v0p1, matrM + 16);
	Mtoeplitz_vm_16x16(rop + 16, vect, m0m1);
	Mtoeplitz_vm_16x16(rop, vect + 16, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
	{
		rop[i] += t0[i];
		rop[i + 16] += t0[i];
	}
}

static inline void M1toeplitz_vm_4x4(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 4
	for(int i = 0; i < 4; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 4
		for(int j = 0; j < 4; j++)
			rop[i] += vect[j] * matr[4 - 1 - j + i];
	}
}

static inline void M1toeplitz_vm_8x8(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	int64_t t0[4], t1[4], t2[4];
	int64_t v0p1[4], m0m1[7], m0m2[7];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 4; i++)
		v0p1[i] = vect[i] + vect[i + 4];
	#pragma GCC unroll 96
	for(int i = 0; i < 7; i++)
	{
		m0m1[i] = matr[i + 4] - matr[i + 8];
		m0m2[i] = matr[i + 4] - matr[i];
	}
	
	M1toeplitz_vm_4x4(t0, v0p1, matr + 4);
	M1toeplitz_vm_4x4(t1, vect, m0m1);
	M1toeplitz_vm_4x4(t2, vect + 4, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 4; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 4] = t0[i] - t1[i];
	}
}

static inline void M1toeplitz_vm_16x16(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	int64_t t0[8], t1[8], t2[8];
	int64_t v0p1[8], m0m1[15], m0m2[15];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
		v0p1[i] = vect[i] + vect[i + 8];
	#pragma GCC unroll 96
	for(int i = 0; i < 15; i++)
	{
		m0m1[i] = matr[i + 8] - matr[i + 16];
		m0m2[i] = matr[i + 8] - matr[i];
	}
	
	M1toeplitz_vm_8x8(t0, v0p1, matr + 8);
	M1toeplitz_vm_8x8(t1, vect, m0m1);
	M1toeplitz_vm_8x8(t2, vect + 8, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 8; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 8] = t0[i] - t1[i];
	}
}

static inline void M1toeplitz_vm_32x32(int64_t* restrict rop, const __int128* restrict vect)
{
	int64_t t0[16], t1[16], t2[16];
	int64_t v0[16], v1[16], v0p1[16], m0m1[31], m0m2[31];
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
	{
		v0p1[i] = vect[i] + vect[i + 16];
		v0[i] = vect[i];
		v1[i] = vect[i + 16];
	}
	#pragma GCC unroll 96
	for(int i = 0; i < 31; i++)
	{
		m0m1[i] = matrM1[i + 16] - matrM1[i + 32];
		m0m2[i] = matrM1[i + 16] - matrM1[i];
	}
	
	M1toeplitz_vm_16x16(t0, v0p1, matrM1 + 16);
	M1toeplitz_vm_16x16(t1, v0, m0m1);
	M1toeplitz_vm_16x16(t2, v1, m0m2);
	
	#pragma GCC unroll 96
	for(int i = 0; i < 16; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 16] = t0[i] - t1[i];
	}
}

static inline void pmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	int64_t matrix[2*N - 1];
	
	_PRAGMAGCCUNROLLLOOP_
	for(int i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];
	
	toeplitz_vm_32x32(R, A->t, matrix);
}

static inline void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Polynomial Montgomery-like method.
	int64_t T[N];
	register uint16_t i;
	
	M1toeplitz_vm_32x32(T, R);
	
	Mtoeplitz_vm_32x32(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

void pmns1664_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the polynomial reduction in a
	// pmns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	pmns_mod_mult_ext_red(R, A, B);
	
	pmns_montg_int_red(res, R);
}


#include <stdint.h>

#include "pmns2048.h"
#include "params2048.h"

const mpnum __P2048__ = &__P__;
const mpnum Gamma2048 = Gi;
const int8_t N2048 = N;
const int8_t RHO2048 = RHO;

static void toeplitz_vm_5x5(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 5
	for(int i = 0; i < 5; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 5
		for(int j = 0; j < 5; j++)
			rop[i] += (__int128) vect[j] * matr[5 - 1 - j + i];
	}
}

static void toeplitz_vm_10x10(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[5], t1[5], t2[5];
	int64_t v0p1[5], m0m1[9], m0m2[9];
	register uint16_t i;
	
	#pragma GCC unroll 5
	for(i = 0; i < 5; i++)
		v0p1[i] = vect[i] + vect[i + 5];
	#pragma GCC unroll 9
	for(i = 0; i < 9; i++)
	{
		m0m1[i] = matr[i + 5] - matr[i + 10];
		m0m2[i] = matr[i + 5] - matr[i];
	}
	
	toeplitz_vm_5x5(t0, v0p1, matr + 5);
	toeplitz_vm_5x5(t1, vect, m0m1);
	toeplitz_vm_5x5(t2, vect + 5, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 5; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 5] = t0[i] - t1[i];
	}
}

static void toeplitz_vm_20x20(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[10], t1[10], t2[10];
	int64_t v0p1[10], m0m1[19], m0m2[19];
	register uint16_t i;
	
	#pragma GCC unroll 10
	for(i = 0; i < 10; i++)
		v0p1[i] = vect[i] + vect[i + 10];
	#pragma GCC unroll 19
	for(i = 0; i < 19; i++)
	{
		m0m1[i] = matr[i + 10] - matr[i + 20];
		m0m2[i] = matr[i + 10] - matr[i];
	}
	
	toeplitz_vm_10x10(t0, v0p1, matr + 10);
	toeplitz_vm_10x10(t1, vect, m0m1);
	toeplitz_vm_10x10(t2, vect + 10, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 10; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 10] = t0[i] - t1[i];
	}
}

static void toeplitz_vm_40x40(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	__int128 t0[20], t1[20], t2[20];
	int64_t v0p1[20], m0m1[39], m0m2[39];
	register uint16_t i;
	
	#pragma GCC unroll 20
	for(i = 0; i < 20; i++)
		v0p1[i] = vect[i] + vect[i + 20];
	#pragma GCC unroll 39
	for(i = 0; i < 39; i++)
	{
		m0m1[i] = matr[i + 20] - matr[i + 40];
		m0m2[i] = matr[i + 20] - matr[i];
	}
	
	toeplitz_vm_20x20(t0, v0p1, matr + 20);
	toeplitz_vm_20x20(t1, vect, m0m1);
	toeplitz_vm_20x20(t2, vect + 20, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 20; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 20] = t0[i] - t1[i];
	}
}

static inline void Movtoeplitz_vm_5x5(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 5
	for(int i = 0; i < 5; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 5
		for(int j = 0; j < 5; j++)
			rop[i] += (__int128) vect[j] * matr[5 - 1 - j + i];
	}
}

static inline void Movtoeplitz_vm_10x10(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	__int128 v0p1[5], t0[5], t1[5], t2[5];
	int64_t m0m1[9], m0m2[9];
	register uint16_t i;
	
	#pragma GCC unroll 5
	for(i = 0; i < 5; i++)
		v0p1[i] = (__int128) vect[i] + vect[i + 5];
	#pragma GCC unroll 9
	for(i = 0; i < 9; i++)
	{
		m0m1[i] = matr[i + 5] - matr[i + 10];
		m0m2[i] = matr[i + 5] - matr[i];
	}
	
	Movtoeplitz_vm_5x5(t0, v0p1, matr + 5);
	Movtoeplitz_vm_5x5(t1, vect, m0m1);
	Movtoeplitz_vm_5x5(t2, vect + 5, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 5; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 5] = t0[i] - t1[i];
	}
}

static void Movtoeplitz_vm_20x20(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
{
	__int128 v0p1[10], t0[10], t1[10], t2[10];
	int64_t m0m1[19], m0m2[19];
	register uint16_t i;
	
	#pragma GCC unroll 10
	for(i = 0; i < 10; i++)
		v0p1[i] = (__int128) vect[i] + vect[i + 10];
	#pragma GCC unroll 19
	for(i = 0; i < 19; i++)
	{
		m0m1[i] = matr[i + 10] - matr[i + 20];
		m0m2[i] = matr[i + 10] - matr[i];
	}
	
	Movtoeplitz_vm_10x10(t0, v0p1, matr + 10);
	Movtoeplitz_vm_10x10(t1, vect, m0m1);
	Movtoeplitz_vm_10x10(t2, vect + 10, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 10; i++)
	{
		rop[i] += t0[i] - t2[i];
		rop[i + 10] += t0[i] - t1[i];
	}
}

static inline void Mtoeplitz_vm_20x20(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 20
	for(int i = 0; i < 20; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 20
		for(int j = 0; j < 20; j++)
			rop[i] += (__int128) vect[j] * matr[20 - 1 - j + i];
	}
}

static void Mtoeplitz_vm_40x40(__int128* restrict rop, const int64_t* restrict vect)
{
	__int128 t0[20], t1[20], t2[20], v0p1[20];
	int64_t m0m1[39], m0m2[39];
	register uint16_t i;
	
	#pragma GCC unroll 20
	for(i = 0; i < 20; i++)
	{
		v0p1[i] = (__int128) vect[i] + vect[i + 20];
	}
	#pragma GCC unroll 39
	for(i = 0; i < 39; i++)
	{
		m0m1[i] = matrM[i + 20] - matrM[i + 40];
		m0m2[i] = matrM[i + 20] - matrM[i];
	}
	
	Movtoeplitz_vm_20x20(t0, v0p1, matrM + 20);
	Mtoeplitz_vm_20x20(t1, vect, m0m1);
	Mtoeplitz_vm_20x20(t2, vect + 20, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 20; i++)
	{
		rop[i] += t0[i] - t2[i];
		rop[i + 20] += t0[i] - t1[i];
	}
}

static inline void M1toeplitz_vm_5x5(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	#pragma GCC unroll 5
	for(int i = 0; i < 5; i++)
	{
		rop[i] = 0;
		#pragma GCC unroll 5
		for(int j = 0; j < 5; j++)
			rop[i] += (__int128) vect[j] * matr[5 - 1 - j + i];
	}
}

static inline void M1toeplitz_vm_10x10(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	int64_t t0[5], t1[5], t2[5];
	int64_t v0p1[5], m0m1[9], m0m2[9];
	register uint16_t i;
	
	#pragma GCC unroll 5
	for(i = 0; i < 5; i++)
	{
		v0p1[i] = vect[i] + vect[i + 5];
	}
	#pragma GCC unroll 9
	for(i = 0; i < 9; i++)
	{
		m0m1[i] = matr[i + 5] - matr[i + 10];
		m0m2[i] = matr[i + 5] - matr[i];
	}
	
	M1toeplitz_vm_5x5(t0, v0p1, matr + 5);
	M1toeplitz_vm_5x5(t1, vect, m0m1);
	M1toeplitz_vm_5x5(t2, vect + 5, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 5; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 5] = t0[i] - t1[i];
	}
}

static inline void M1toeplitz_vm_20x20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	int64_t t0[10], t1[10], t2[10];
	int64_t v0p1[10], m0m1[19], m0m2[19];
	register uint16_t i;
	
	#pragma GCC unroll 10
	for(i = 0; i < 10; i++)
	{
		v0p1[i] = vect[i] + vect[i + 10];
	}
	#pragma GCC unroll 19
	for(i = 0; i < 19; i++)
	{
		m0m1[i] = matr[i + 10] - matr[i + 20];
		m0m2[i] = matr[i + 10] - matr[i];
	}
	
	M1toeplitz_vm_10x10(t0, v0p1, matr + 10);
	M1toeplitz_vm_10x10(t1, vect, m0m1);
	M1toeplitz_vm_10x10(t2, vect + 10, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 10; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 10] = t0[i] - t1[i];
	}
}

static void M1toeplitz_vm_40x40(int64_t* restrict rop, const __int128* restrict vect)
{
	int64_t t0[20], t1[20], t2[20];
	int64_t v0p1[20], m0m1[39], m0m2[39], v0[20], v1[20];
	register uint16_t i;
	
	#pragma GCC unroll 20
	for(i = 0; i < 20; i++)
	{
		v0[i] = vect[i];
		v1[i] = vect[i + 20];
		v0p1[i] = vect[i] + vect[i + 20];
	}
	#pragma GCC unroll 39
	for(i = 0; i < 39; i++)
	{
		m0m1[i] = matrM1[i + 20] - matrM1[i + 40];
		m0m2[i] = matrM1[i + 20] - matrM1[i];
	}
	
	M1toeplitz_vm_20x20(t0, v0p1, matrM1 + 20);
	M1toeplitz_vm_20x20(t1, v0, m0m1);
	M1toeplitz_vm_20x20(t2, v1, m0m2);
	
	#pragma GCC unroll 96
	for(i = 0; i < 20; i++)
	{
		rop[i] = t0[i] - t2[i];
		rop[i + 20] = t0[i] - t1[i];
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
	
	toeplitz_vm_40x40(R, A->t, matrix);
}

static inline void pmns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Polynomial Montgomery-like method.
	int64_t T[N];
	register uint16_t i;
	
	M1toeplitz_vm_40x40(T, R);
	
	Mtoeplitz_vm_40x40(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

void pmns2048_montg_mult(restrict poly res, const restrict poly A,
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


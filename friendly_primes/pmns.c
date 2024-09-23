#include "pmns.h"

#define SCHOOLBOOK(X) {\
	for(int i = 0; i < X; i++)\
	{\
		rop[i] = 0;\
		for(int j = 0; j < X; j++)\
			rop[i] += (__int128) vect[j] * matr[X - 1 - j + i];\
	}\
}

#define TOEP22TOP(X, F) {\
	__int128 t0[X/2], t1[X/2], t2[X/2];\
	int64_t v0p1[X/2], m0m1[X-1], m0m2[X-1];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = vect[i] + vect[i + X/2];\
	}\
	for(int i = 0; i < X-1; i++)\
	{\
		m0m1[i] = matr[i + X/2] - matr[i + X];\
		m0m2[i] = matr[i + X/2] - matr[i];\
	}\
	F (t0, v0p1, matr + X/2);\
	F (t1, vect, m0m1);\
	F (t2, vect + X/2, m0m2);\
	for(int i = 0; i < X/2; i++)\
	{\
		rop[i] = t0[i] - t2[i];\
		rop[i + X/2] = t0[i] - t1[i];\
	}\
}

#define TOEP33TOP(X, F) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3];\
	int64_t m03, v1m2[X/3], v0m2[X/3], v0m1[X/3], m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = vect[i] - vect[X/3 + i];\
		v0m2[i] = vect[i] - vect[2*X/3 + i];\
	}\
	for(int i = 0; i < (2*X/3) - 1; i++)\
	{\
		m03 = matr[i + 2*X/3] + matr[i + X/3];\
		m034[i] = m03 + matr[i];\
		m013[i] = m03 +  matr[i + X];\
		m012[i] = matr[i + 2*X/3] +  matr[i + X] + matr[i + 4*X/3];\
	}\
	F (t0, vect + 2*X/3, m034);\
	F (t1, vect + X/3, m013);\
	F (t2, vect, m012);\
	F (t3, v1m2,  matr + X/3);\
	F (t4, v0m2, matr + 2*X/3);\
	F (t5, v0m1,  matr + X);\
	for(int i = 0; i < X/3; i++)\
	{\
		rop[i] = t0[i] + t3[i] + t4[i];\
		rop[i + X/3] = t1[i] - t3[i] + t5[i];\
		rop[i + 2*X/3] = t2[i] - t4[i] - t5[i];\
	}\
}

#define pmns_toeplitz_ext_red(N) \
{\
	int64_t matr[2*N - 1];\
	for(int i = 0; i < N-1; i++)\
	{\
		matr[i + N - 1] = B->t[i];\
		matr[i] = B->t[1 + i] * LAMBDA;\
	}\
	matr[2*N - 2] = B->t[N - 1];\
	toeplitz_vm(R, A->t, matr);\
}

#define pmns_normal_ext_red(N) \
{\
	int64_t btab[N];\
	for(int i = 1; i < N; i++)\
		btab[N-i] = B->t[i]*LAMBDA;\
	for(int i = 0; i < N; i++)\
	{\
		R[i] = (__int128) A->t[0] * B->t[i];\
		for(int j = 1; j < i + 1; j++)\
			R[i] += (__int128) A->t[j] * B->t[i - j];\
	}\
	for(int i = 0; i < N - 1; i++)\
	{\
		for(int j = 1; j < N - i; j++)\
			R[i] += (__int128) A->t[i + j] * btab[j];\
	}\
}

#define pmns_mult_by_g(N) \
{\
	int64_t *pT = (int64_t*)T;\
	R[0] += (__int128)pT[0] * GAMMA - (__int128)pT[N - 1] * LAMED;\
	for(int i = 1; i < N - 1; i++)\
	{\
		R[i] += (__int128)pT[i] * GAMMA - pT[i - 1];\
	}\
	R[N - 1] += (__int128)pT[N - 1] * GIMEL - pT[N - 2];\
}

#define pmns_mult_by_linearg1(N) \
{\
	uint64_t Z = 0;\
	for(int j = 0; j < N; j++)\
		Z += (uint64_t)R[j] * (uint64_t)lastcol[j];\
	Z = -Z;\
	T[N - 1] = Z;\
	for(int i = 0; i < N - 1; i++)\
	{\
		Z *= GAMMA;\
		Z += (uint64_t)R[N - 1 - i];\
		T[N - 2 - i] = Z;\
	}\
}

#define pmns_mult_by_sparseg1(N) \
{\
	for(int i = 0; i < N - 2; i++)\
	{\
		T[i] = (uint64_t)R[i + 1] + (uint64_t)R[i + 2] * GAMMA;\
	}\
	T[N - 2] = (uint64_t)R[0] * GAMMALAMM1 + (uint64_t)R[N - 1];\
	T[N - 1] = (uint64_t)R[0] * ONELAMM1 + (uint64_t)R[1] * GAMMALAMM2;\
}

#define pmns_linear_int_red(N) \
{\
	uint64_t T[N];\
	pmns_mult_by_linearg1(N)\
	pmns_mult_by_g(N)\
	for(int i = 0; i < N; i++)\
		res->t[i] = (int64_t)(R[i] >> 64);\
}

#define pmns_sparse_int_red(N) \
{\
	uint64_t T[N];\
	pmns_mult_by_sparseg1(N)\
	pmns_mult_by_g(N)\
	for(int i = 0; i < N; i++)\
		res->t[i] = (int64_t)(R[i] >> 64);\
}

#define _pmns_linearred_mult(N) \
{\
	__int128 R[N];\
	pmns_normal_ext_red(N)\
	pmns_linear_int_red(N)\
}

#define _pmns_doublesparse_mult(N) \
{\
	__int128 R[N];\
	pmns_normal_ext_red(N)\
	pmns_sparse_int_red(N)\
}

#define _pmns_doublesparse_toep(N) \
{\
	__int128 R[N];\
	pmns_toeplitz_ext_red(N)\
	pmns_sparse_int_red(N)\
}

#define _pmns_linearred_toep(N) \
{\
	__int128 R[N];\
	pmns_toeplitz_ext_red(N)\
	pmns_linear_int_red(N)\
}

uint16_t N;
uint8_t ALPHA = 1;
uint64_t GAMMA, GIMEL;
int16_t LAMBDA, LAMED;
uint64_t lastcol[144];
uint64_t GAMMALAMM1, GAMMALAMM2, ONELAMM1;

static inline void schoolbook4(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(4)

static inline void schoolbook3(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(3)

#ifdef VARIOUSNVALUES

static inline void schoolbook5(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(5)

static inline void schoolbook6(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(6)

static inline void schoolbook7(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(7)

#endif

#ifdef LARGEPMNS

static inline void toeplitz9(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP33TOP(9, schoolbook3)

static inline void toeplitz18(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP22TOP(18, toeplitz9)

static inline void toeplitz36(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP22TOP(36, toeplitz18)

static inline void toeplitz72(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP22TOP(72, toeplitz36)

#endif

static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
{
	if(N == 8)
		TOEP22TOP(8, schoolbook4)
	else if(N == 9)
		TOEP33TOP(9, schoolbook3)
	#ifdef VARIOUSNVALUES
	else if(N == 10)
		TOEP22TOP(10, schoolbook5)
	else if(N == 12)
		TOEP22TOP(12, schoolbook6)
	else if(N == 14)
		TOEP22TOP(14, schoolbook7)
	#endif
	#ifdef LARGEPMNS
	else if(N == 18)
		toeplitz18(rop, vect, matr);
	else if(N == 36)
		toeplitz36(rop, vect, matr);
	else if(N == 72)
		toeplitz72(rop, vect, matr);
	else if(N == 144)
		TOEP22TOP(144, toeplitz72)
	#endif
}

void setpmnsalpha(uint8_t newALPHA)
{
	ALPHA = newALPHA;
}

void setpmnsparams(uint16_t newN, uint64_t newGAMMA, int16_t newLAMBDA)
{
	N = newN;
	GAMMA = newGAMMA;
	GIMEL = newGAMMA;
	LAMBDA = newLAMBDA;
	LAMED = newLAMBDA;
}

void setpmnspsi(uint64_t psi)
{
	GIMEL /= psi;
	LAMED /= psi;
}

void setpmnsaux(uint64_t newGIMEL, int16_t newLAMED)
{
	GIMEL = newGIMEL;
	LAMED = newLAMED;
}

void setpmnslastcol(uint64_t newcol[])
{
	for(int i = 0; i < N; i++)
		lastcol[i] = newcol[i];
}

void setpmnssparseinv(uint64_t newGL1, uint64_t newGL2, uint64_t newOL1)
{
	GAMMALAMM1 = newGL1;
	GAMMALAMM2 = newGL2;
	ONELAMM1 = newOL1;
}

static inline void _pmns_linearred_nine(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t a0, a1, a2, a3, a4, a5, a6, a7, a8, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;
	__int128 r0, r1, r2, r3, r4, r5, r6, r7, r8;
	int64_t v1m2_0, v1m2_1, v1m2_2, v0m2_0, v0m2_1, v0m2_2, v0m1_0, v0m1_1, v0m1_2, m034_0, m034_1, m034_2, m034_3, m034_4, m013_0, m013_1, m013_2, m013_3, m013_4, m012_0, m012_1, m012_2, m012_3, m012_4;
	__int128 t3_0, t3_1, t3_2, t4_0, t4_1, t4_2, t5_0, t5_1, t5_2;
	int64_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
	
	a0 = A->t[0];
	a1 = A->t[1];
	a2 = A->t[2];
	a3 = A->t[3];
	a4 = A->t[4];
	a5 = A->t[5];
	a6 = A->t[6];
	a7 = A->t[7];
	a8 = A->t[8];
	
	b0 = B->t[1] * LAMBDA;
	b1 = B->t[2] * LAMBDA;
	b2 = B->t[3] * LAMBDA;
	b3 = B->t[4] * LAMBDA;
	b4 = B->t[5] * LAMBDA;
	b5 = B->t[6] * LAMBDA;
	b6 = B->t[7] * LAMBDA;
	b7 = B->t[8] * LAMBDA;
	b8 = B->t[0];
	b9 = B->t[1];
	b10 = B->t[2];
	b11 = B->t[3];
	b12 = B->t[4];
	b13 = B->t[5];
	b14 = B->t[6];
	b15 = B->t[7];
	b16 = B->t[8];
	
	v1m2_0 = a3 - a6;
	v1m2_1 = a4 - a7;
	v1m2_2 = a5 - a8;
	v0m1_0 = a0 - a3;
	v0m1_1 = a1 - a4;
	v0m1_2 = a2 - a5;
	v0m2_0 = a0 - a6;
	v0m2_1 = a1 - a7;
	v0m2_2 = a2 - a8;
	
	m034_0 = b6 + b3 + b0;
	m013_0 = b6 + b3 + b9;
	m012_0 = b6 + b9 + b12;
	m034_1 = b7 + b4 + b1;
	m013_1 = b7 + b4 + b10;
	m012_1 = b7 + b10 + b13;
	m034_2 = b8 + b5 + b2;
	m013_2 = b8 + b5 + b11;
	m012_2 = b8 + b11 + b14;
	m034_3 = b9 + b6 + b3;
	m013_3 = b9 + b6 + b12;
	m012_3 = b9 + b12 + b15;
	m034_4 = b10 + b7 + b4;
	m013_4 = b10 + b7 + b13;
	m012_4 = b10 + b13 + b16;
	
	t3_0 = (__int128) v1m2_0 * b5 + (__int128) v1m2_1 * b4 + (__int128) v1m2_2 * b3;
	t3_1 = (__int128) v1m2_0 * b6 + (__int128) v1m2_1 * b5 + (__int128) v1m2_2 * b4;
	t3_2 = (__int128) v1m2_0 * b7 + (__int128) v1m2_1 * b6 + (__int128) v1m2_2 * b5;
	t4_0 = (__int128) v0m2_0 * b8 + (__int128) v0m2_1 * b7 + (__int128) v0m2_2 * b6;
	t4_1 = (__int128) v0m2_0 * b9 + (__int128) v0m2_1 * b8 + (__int128) v0m2_2 * b7;
	t4_2 = (__int128) v0m2_0 * b10 + (__int128) v0m2_1 * b9 + (__int128) v0m2_2 * b8;
	t5_0 = (__int128) v0m1_0 * b11 + (__int128) v0m1_1 * b10 + (__int128) v0m1_2 * b9;
	t5_1 = (__int128) v0m1_0 * b12 + (__int128) v0m1_1 * b11 + (__int128) v0m1_2 * b10;
	t5_2 = (__int128) v0m1_0 * b13 + (__int128) v0m1_1 * b12 + (__int128) v0m1_2 * b11;
	
	r0 = (__int128) a6 * m034_2 + (__int128) a7 * m034_1 + (__int128) a8 * m034_0 + t3_0 + t4_0;
	r1 = (__int128) a6 * m034_3 + (__int128) a7 * m034_2 + (__int128) a8 * m034_1 + t3_1 + t4_1;
	r2 = (__int128) a6 * m034_4 + (__int128) a7 * m034_3 + (__int128) a8 * m034_2 + t3_2 + t4_2;
	r3 = (__int128) a3 * m013_2 + (__int128) a4 * m013_1 + (__int128) a5 * m013_0 + t5_0 - t3_0;
	r4 = (__int128) a3 * m013_3 + (__int128) a4 * m013_2 + (__int128) a5 * m013_1 + t5_1 - t3_1;
	r5 = (__int128) a3 * m013_4 + (__int128) a4 * m013_3 + (__int128) a5 * m013_2 + t5_2 - t3_2;
	r6 = (__int128) a0 * m012_2 + (__int128) a1 * m012_1 + (__int128) a2 * m012_0 - t4_0 - t5_0;
	r7 = (__int128) a0 * m012_3 + (__int128) a1 * m012_2 + (__int128) a2 * m012_1 - t4_1 - t5_1;
	r8 = (__int128) a0 * m012_4 + (__int128) a1 * m012_3 + (__int128) a2 * m012_2 - t4_2 - t5_2;
	
	uint64_t Z = (uint64_t)r0 * (uint64_t) lastcol[0]
							+(uint64_t)r1 * (uint64_t) lastcol[1]
							+(uint64_t)r2 * (uint64_t) lastcol[2]
							+(uint64_t)r3 * (uint64_t) lastcol[3]
							+(uint64_t)r4 * (uint64_t) lastcol[4]
							+(uint64_t)r5 * (uint64_t) lastcol[5]
							+(uint64_t)r6 * (uint64_t) lastcol[6]
							+(uint64_t)r7 * (uint64_t) lastcol[7]
							+(uint64_t)r8 * (uint64_t) lastcol[8];
	
	t8 = Z;
	Z = Z * GAMMA - (uint64_t) r8; t7 = Z;
	Z = Z * GAMMA - (uint64_t) r7; t6 = Z;
	Z = Z * GAMMA - (uint64_t) r6; t5 = Z;
	Z = Z * GAMMA - (uint64_t) r5; t4 = Z;
	Z = Z * GAMMA - (uint64_t) r4; t3 = Z;
	Z = Z * GAMMA - (uint64_t) r3; t2 = Z;
	Z = Z * GAMMA - (uint64_t) r2; t1 = Z;
	Z = Z * GAMMA - (uint64_t) r1; t0 = Z;
	
	r0 += (__int128)t8 * LAMED - (__int128) t0 * GAMMA;
	r1 += t0 - (__int128) t1 * GAMMA;
	r2 += t1 - (__int128) t2 * GAMMA;
	r3 += t2 - (__int128) t3 * GAMMA;
	r4 += t3 - (__int128) t4 * GAMMA;
	r5 += t4 - (__int128) t5 * GAMMA;
	r6 += t5 - (__int128) t6 * GAMMA;
	r7 += t6 - (__int128) t7 * GAMMA;
	r8 += t7 - (__int128) t8 * GIMEL;
	
	res->t[0] = (int64_t)(r0>>64);
	res->t[1] = (int64_t)(r1>>64);
	res->t[2] = (int64_t)(r2>>64);
	res->t[3] = (int64_t)(r3>>64);
	res->t[4] = (int64_t)(r4>>64);
	res->t[5] = (int64_t)(r5>>64);
	res->t[6] = (int64_t)(r6>>64);
	res->t[7] = (int64_t)(r7>>64);
	res->t[8] = (int64_t)(r8>>64);
}

static inline void _pmns_bino_linear_nine(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t a0, a1, a2, a3, a4, a5, a6, a7, a8, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;
	__int128 r0, r1, r2, r3, r4, r5, r6, r7, r8;
	int64_t v1m2_0, v1m2_1, v1m2_2, v0m2_0, v0m2_1, v0m2_2, v0m1_0, v0m1_1, v0m1_2, m034_0, m034_1, m034_2, m034_3, m034_4, m013_0, m013_1, m013_2, m013_3, m013_4, m012_0, m012_1, m012_2, m012_3, m012_4;
	__int128 t3_0, t3_1, t3_2, t4_0, t4_1, t4_2, t5_0, t5_1, t5_2;
	int64_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
	
	a0 = A->t[0];
	a1 = A->t[1];
	a2 = A->t[2];
	a3 = A->t[3];
	a4 = A->t[4];
	a5 = A->t[5];
	a6 = A->t[6];
	a7 = A->t[7];
	a8 = A->t[8];
	
	b0 = B->t[1] * LAMBDA;
	b1 = B->t[2] * LAMBDA;
	b2 = B->t[3] * LAMBDA;
	b3 = B->t[4] * LAMBDA;
	b4 = B->t[5] * LAMBDA;
	b5 = B->t[6] * LAMBDA;
	b6 = B->t[7] * LAMBDA;
	b7 = B->t[8] * LAMBDA;
	b8 = B->t[0] * ALPHA;
	b9 = B->t[1] * ALPHA;
	b10 = B->t[2] * ALPHA;
	b11 = B->t[3] * ALPHA;
	b12 = B->t[4] * ALPHA;
	b13 = B->t[5] * ALPHA;
	b14 = B->t[6] * ALPHA;
	b15 = B->t[7] * ALPHA;
	b16 = B->t[8] * ALPHA;
	
	v1m2_0 = a3 - a6;
	v1m2_1 = a4 - a7;
	v1m2_2 = a5 - a8;
	v0m1_0 = a0 - a3;
	v0m1_1 = a1 - a4;
	v0m1_2 = a2 - a5;
	v0m2_0 = a0 - a6;
	v0m2_1 = a1 - a7;
	v0m2_2 = a2 - a8;
	
	m034_0 = b6 + b3 + b0;
	m013_0 = b6 + b3 + b9;
	m012_0 = b6 + b9 + b12;
	m034_1 = b7 + b4 + b1;
	m013_1 = b7 + b4 + b10;
	m012_1 = b7 + b10 + b13;
	m034_2 = b8 + b5 + b2;
	m013_2 = b8 + b5 + b11;
	m012_2 = b8 + b11 + b14;
	m034_3 = b9 + b6 + b3;
	m013_3 = b9 + b6 + b12;
	m012_3 = b9 + b12 + b15;
	m034_4 = b10 + b7 + b4;
	m013_4 = b10 + b7 + b13;
	m012_4 = b10 + b13 + b16;
	
	t3_0 = (__int128) v1m2_0 * b5 + (__int128) v1m2_1 * b4 + (__int128) v1m2_2 * b3;
	t3_1 = (__int128) v1m2_0 * b6 + (__int128) v1m2_1 * b5 + (__int128) v1m2_2 * b4;
	t3_2 = (__int128) v1m2_0 * b7 + (__int128) v1m2_1 * b6 + (__int128) v1m2_2 * b5;
	t4_0 = (__int128) v0m2_0 * b8 + (__int128) v0m2_1 * b7 + (__int128) v0m2_2 * b6;
	t4_1 = (__int128) v0m2_0 * b9 + (__int128) v0m2_1 * b8 + (__int128) v0m2_2 * b7;
	t4_2 = (__int128) v0m2_0 * b10 + (__int128) v0m2_1 * b9 + (__int128) v0m2_2 * b8;
	t5_0 = (__int128) v0m1_0 * b11 + (__int128) v0m1_1 * b10 + (__int128) v0m1_2 * b9;
	t5_1 = (__int128) v0m1_0 * b12 + (__int128) v0m1_1 * b11 + (__int128) v0m1_2 * b10;
	t5_2 = (__int128) v0m1_0 * b13 + (__int128) v0m1_1 * b12 + (__int128) v0m1_2 * b11;
	
	r0 = (__int128) a6 * m034_2 + (__int128) a7 * m034_1 + (__int128) a8 * m034_0 + t3_0 + t4_0;
	r1 = (__int128) a6 * m034_3 + (__int128) a7 * m034_2 + (__int128) a8 * m034_1 + t3_1 + t4_1;
	r2 = (__int128) a6 * m034_4 + (__int128) a7 * m034_3 + (__int128) a8 * m034_2 + t3_2 + t4_2;
	r3 = (__int128) a3 * m013_2 + (__int128) a4 * m013_1 + (__int128) a5 * m013_0 + t5_0 - t3_0;
	r4 = (__int128) a3 * m013_3 + (__int128) a4 * m013_2 + (__int128) a5 * m013_1 + t5_1 - t3_1;
	r5 = (__int128) a3 * m013_4 + (__int128) a4 * m013_3 + (__int128) a5 * m013_2 + t5_2 - t3_2;
	r6 = (__int128) a0 * m012_2 + (__int128) a1 * m012_1 + (__int128) a2 * m012_0 - t4_0 - t5_0;
	r7 = (__int128) a0 * m012_3 + (__int128) a1 * m012_2 + (__int128) a2 * m012_1 - t4_1 - t5_1;
	r8 = (__int128) a0 * m012_4 + (__int128) a1 * m012_3 + (__int128) a2 * m012_2 - t4_2 - t5_2;
	
	uint64_t Z = (uint64_t)r0 * (uint64_t) lastcol[0]
							+(uint64_t)r1 * (uint64_t) lastcol[1]
							+(uint64_t)r2 * (uint64_t) lastcol[2]
							+(uint64_t)r3 * (uint64_t) lastcol[3]
							+(uint64_t)r4 * (uint64_t) lastcol[4]
							+(uint64_t)r5 * (uint64_t) lastcol[5]
							+(uint64_t)r6 * (uint64_t) lastcol[6]
							+(uint64_t)r7 * (uint64_t) lastcol[7]
							+(uint64_t)r8 * (uint64_t) lastcol[8];
	
	t8 = Z;
	Z = Z * ALPHA * GAMMA - (uint64_t) r8; t7 = Z;
	Z = Z * GAMMA - (uint64_t) r7; t6 = Z;
	Z = Z * GAMMA - (uint64_t) r6; t5 = Z;
	Z = Z * GAMMA - (uint64_t) r5; t4 = Z;
	Z = Z * GAMMA - (uint64_t) r4; t3 = Z;
	Z = Z * GAMMA - (uint64_t) r3; t2 = Z;
	Z = Z * GAMMA - (uint64_t) r2; t1 = Z;
	Z = Z * GAMMA - (uint64_t) r1; t0 = Z;
	
	r0 += (__int128)t8 * LAMED - (__int128) t0 * GAMMA;
	r1 += t0 - (__int128) t1 * GAMMA;
	r2 += t1 - (__int128) t2 * GAMMA;
	r3 += t2 - (__int128) t3 * GAMMA;
	r4 += t3 - (__int128) t4 * GAMMA;
	r5 += t4 - (__int128) t5 * GAMMA;
	r6 += t5 - (__int128) t6 * GAMMA;
	r7 += t6 - (__int128) t7 * GAMMA;
	r8 += t7 - (__int128) t8 * GIMEL;
	
	res->t[0] = (int64_t)(r0>>64);
	res->t[1] = (int64_t)(r1>>64);
	res->t[2] = (int64_t)(r2>>64);
	res->t[3] = (int64_t)(r3>>64);
	res->t[4] = (int64_t)(r4>>64);
	res->t[5] = (int64_t)(r5>>64);
	res->t[6] = (int64_t)(r6>>64);
	res->t[7] = (int64_t)(r7>>64);
	res->t[8] = (int64_t)(r8>>64);
}

static inline void _pmns_doublesparse_five(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t a0, a1, a2, a3, a4, b1, b2, b3, b4, b5, b6, b7, b8, b9;
	__int128 r0, r1, r2, r3, r4;
	int64_t v0p1_0, v0p1_1, m0m1_0, m0m1_1, m0m1_2, m0m1_3, m0m2_1, m0m2_2, m0m2_3, m0m2_4;
	__int128 t0_0, t0_1, t0_2;
	int64_t t0, t1, t2, t3, t4;
	
	a0 = A->t[0];
	a1 = A->t[1];
	a2 = A->t[2];
	a3 = A->t[3];
	a4 = A->t[4];
	
	b1 = B->t[1] * LAMBDA;
	b2 = B->t[2] * LAMBDA;
	b3 = B->t[3] * LAMBDA;
	b4 = B->t[4] * LAMBDA;
	b5 = B->t[0];
	b6 = B->t[1];
	b7 = B->t[2];
	b8 = B->t[3];
	b9 = B->t[4];
	
	v0p1_0 = a0 + a3;
	v0p1_1 = a1 + a4;
	
	m0m1_0 = b3 - b6;
	m0m1_1 = b4 - b7;
	m0m1_2 = b5 - b8;
	m0m1_3 = b6 - b9;
	m0m2_1 = b4 - b1;
	m0m2_2 = b5 - b2;
	m0m2_3 = b6 - b3;
	m0m2_4 = b7 - b4;
	
	t0_0 = (__int128) v0p1_0 * b5 + (__int128) v0p1_1 * b4 + (__int128) a2 * b3;
	t0_1 = (__int128) v0p1_0 * b6 + (__int128) v0p1_1 * b5 + (__int128) a2 * b4;
	t0_2 = (__int128) v0p1_0 * b7 + (__int128) v0p1_1 * b6 + (__int128) a2 * b5;
	
	r0 = t0_0 - ((__int128) a3 * m0m2_2 + (__int128) a4 * m0m2_1);
	r1 = t0_1 - ((__int128) a3 * m0m2_3 + (__int128) a4 * m0m2_2);
	r2 = t0_2 - ((__int128) a3 * m0m2_4 + (__int128) a4 * m0m2_3);
	r3 = t0_0 - ((__int128) a0 * m0m1_2 + (__int128) a1 * m0m1_1 + (__int128) a2 * m0m1_0);
	r4 = t0_1 - ((__int128) a0 * m0m1_3 + (__int128) a1 * m0m1_2 + (__int128) a2 * m0m1_1);
	
	t0 = ((uint64_t)r1 + (uint64_t)((uint64_t)r2*GAMMA));
	t1 = ((uint64_t)r2 + (uint64_t)((uint64_t)r3*GAMMA));
	t2 = ((uint64_t)r3 + (uint64_t)((uint64_t)r4*GAMMA));
	t3 = ((uint64_t)r4 + (uint64_t)((uint64_t)r0 * GAMMALAMM1));
	t4 = ((uint64_t)((uint64_t)r0 * ONELAMM1) + (uint64_t)((uint64_t)r1 * GAMMALAMM2));
	
	r0 += (__int128) t0 * GAMMA - (__int128)t4 * LAMED;
	r1 += (__int128) t1 * GAMMA - t0;
	r2 += (__int128) t2 * GAMMA - t1;
	r3 += (__int128) t3 * GAMMA - t2;
	r4 += (__int128) t4 * GIMEL - t3;
	
	res->t[0] = (int64_t)(r0>>64);
	res->t[1] = (int64_t)(r1>>64);
	res->t[2] = (int64_t)(r2>>64);
	res->t[3] = (int64_t)(r3>>64);
	res->t[4] = (int64_t)(r4>>64);
}

static inline void _pmns_doublesparse_nine(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t a0, a1, a2, a3, a4, a5, a6, a7, a8, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;
	__int128 r0, r1, r2, r3, r4, r5, r6, r7, r8;
	int64_t v1m2_0, v1m2_1, v1m2_2, v0m2_0, v0m2_1, v0m2_2, v0m1_0, v0m1_1, v0m1_2, m034_0, m034_1, m034_2, m034_3, m034_4, m013_0, m013_1, m013_2, m013_3, m013_4, m012_0, m012_1, m012_2, m012_3, m012_4;
	__int128 t3_0, t3_1, t3_2, t4_0, t4_1, t4_2, t5_0, t5_1, t5_2;
	int64_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
	
	a0 = A->t[0];
	a1 = A->t[1];
	a2 = A->t[2];
	a3 = A->t[3];
	a4 = A->t[4];
	a5 = A->t[5];
	a6 = A->t[6];
	a7 = A->t[7];
	a8 = A->t[8];
	
	b0 = B->t[1] * LAMBDA;
	b1 = B->t[2] * LAMBDA;
	b2 = B->t[3] * LAMBDA;
	b3 = B->t[4] * LAMBDA;
	b4 = B->t[5] * LAMBDA;
	b5 = B->t[6] * LAMBDA;
	b6 = B->t[7] * LAMBDA;
	b7 = B->t[8] * LAMBDA;
	b8 = B->t[0];
	b9 = B->t[1];
	b10 = B->t[2];
	b11 = B->t[3];
	b12 = B->t[4];
	b13 = B->t[5];
	b14 = B->t[6];
	b15 = B->t[7];
	b16 = B->t[8];
	
	v1m2_0 = a3 - a6;
	v1m2_1 = a4 - a7;
	v1m2_2 = a5 - a8;
	v0m1_0 = a0 - a3;
	v0m1_1 = a1 - a4;
	v0m1_2 = a2 - a5;
	v0m2_0 = a0 - a6;
	v0m2_1 = a1 - a7;
	v0m2_2 = a2 - a8;
	
	m034_0 = b6 + b3 + b0;
	m013_0 = b6 + b3 + b9;
	m012_0 = b6 + b9 + b12;
	m034_1 = b7 + b4 + b1;
	m013_1 = b7 + b4 + b10;
	m012_1 = b7 + b10 + b13;
	m034_2 = b8 + b5 + b2;
	m013_2 = b8 + b5 + b11;
	m012_2 = b8 + b11 + b14;
	m034_3 = b9 + b6 + b3;
	m013_3 = b9 + b6 + b12;
	m012_3 = b9 + b12 + b15;
	m034_4 = b10 + b7 + b4;
	m013_4 = b10 + b7 + b13;
	m012_4 = b10 + b13 + b16;
	
	t3_0 = (__int128) v1m2_0 * b5 + (__int128) v1m2_1 * b4 + (__int128) v1m2_2 * b3;
	t3_1 = (__int128) v1m2_0 * b6 + (__int128) v1m2_1 * b5 + (__int128) v1m2_2 * b4;
	t3_2 = (__int128) v1m2_0 * b7 + (__int128) v1m2_1 * b6 + (__int128) v1m2_2 * b5;
	t4_0 = (__int128) v0m2_0 * b8 + (__int128) v0m2_1 * b7 + (__int128) v0m2_2 * b6;
	t4_1 = (__int128) v0m2_0 * b9 + (__int128) v0m2_1 * b8 + (__int128) v0m2_2 * b7;
	t4_2 = (__int128) v0m2_0 * b10 + (__int128) v0m2_1 * b9 + (__int128) v0m2_2 * b8;
	t5_0 = (__int128) v0m1_0 * b11 + (__int128) v0m1_1 * b10 + (__int128) v0m1_2 * b9;
	t5_1 = (__int128) v0m1_0 * b12 + (__int128) v0m1_1 * b11 + (__int128) v0m1_2 * b10;
	t5_2 = (__int128) v0m1_0 * b13 + (__int128) v0m1_1 * b12 + (__int128) v0m1_2 * b11;
	
	r0 = (__int128) a6 * m034_2 + (__int128) a7 * m034_1 + (__int128) a8 * m034_0 + t3_0 + t4_0;
	r1 = (__int128) a6 * m034_3 + (__int128) a7 * m034_2 + (__int128) a8 * m034_1 + t3_1 + t4_1;
	r2 = (__int128) a6 * m034_4 + (__int128) a7 * m034_3 + (__int128) a8 * m034_2 + t3_2 + t4_2;
	r3 = (__int128) a3 * m013_2 + (__int128) a4 * m013_1 + (__int128) a5 * m013_0 + t5_0 - t3_0;
	r4 = (__int128) a3 * m013_3 + (__int128) a4 * m013_2 + (__int128) a5 * m013_1 + t5_1 - t3_1;
	r5 = (__int128) a3 * m013_4 + (__int128) a4 * m013_3 + (__int128) a5 * m013_2 + t5_2 - t3_2;
	r6 = (__int128) a0 * m012_2 + (__int128) a1 * m012_1 + (__int128) a2 * m012_0 - t4_0 - t5_0;
	r7 = (__int128) a0 * m012_3 + (__int128) a1 * m012_2 + (__int128) a2 * m012_1 - t4_1 - t5_1;
	r8 = (__int128) a0 * m012_4 + (__int128) a1 * m012_3 + (__int128) a2 * m012_2 - t4_2 - t5_2;
	
	t0 = ((uint64_t)r1 + (uint64_t)((uint64_t)r2*GAMMA));
	t1 = ((uint64_t)r2 + (uint64_t)((uint64_t)r3*GAMMA));
	t2 = ((uint64_t)r3 + (uint64_t)((uint64_t)r4*GAMMA));
	t3 = ((uint64_t)r4 + (uint64_t)((uint64_t)r5*GAMMA));
	t4 = ((uint64_t)r5 + (uint64_t)((uint64_t)r6*GAMMA));
	t5 = ((uint64_t)r6 + (uint64_t)((uint64_t)r7*GAMMA));
	t6 = ((uint64_t)r7 + (uint64_t)((uint64_t)r8*GAMMA));
	t7 = ((uint64_t)r8 + (uint64_t)((uint64_t)r0 * GAMMALAMM1));
	t8 = ((uint64_t)((uint64_t)r0 * ONELAMM1) + (uint64_t)((uint64_t)r1 * GAMMALAMM2));
	
	r0 += (__int128) t0 * GAMMA - (__int128)t8 * LAMED;
	r1 += (__int128) t1 * GAMMA - t0;
	r2 += (__int128) t2 * GAMMA - t1;
	r3 += (__int128) t3 * GAMMA - t2;
	r4 += (__int128) t4 * GAMMA - t3;
	r5 += (__int128) t5 * GAMMA - t4;
	r6 += (__int128) t6 * GAMMA - t5;
	r7 += (__int128) t7 * GAMMA - t6;
	r8 += (__int128) t8 * GIMEL - t7;
	
	res->t[0] = (int64_t)(r0>>64);
	res->t[1] = (int64_t)(r1>>64);
	res->t[2] = (int64_t)(r2>>64);
	res->t[3] = (int64_t)(r3>>64);
	res->t[4] = (int64_t)(r4>>64);
	res->t[5] = (int64_t)(r5>>64);
	res->t[6] = (int64_t)(r6>>64);
	res->t[7] = (int64_t)(r7>>64);
	res->t[8] = (int64_t)(r8>>64);
}

static inline void _pmns_bino_sparse_nine(restrict poly res, const restrict poly A, const restrict poly B)
{
	int64_t a0, a1, a2, a3, a4, a5, a6, a7, a8, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;
	__int128 r0, r1, r2, r3, r4, r5, r6, r7, r8;
	int64_t v1m2_0, v1m2_1, v1m2_2, v0m2_0, v0m2_1, v0m2_2, v0m1_0, v0m1_1, v0m1_2, m034_0, m034_1, m034_2, m034_3, m034_4, m013_0, m013_1, m013_2, m013_3, m013_4, m012_0, m012_1, m012_2, m012_3, m012_4;
	__int128 t3_0, t3_1, t3_2, t4_0, t4_1, t4_2, t5_0, t5_1, t5_2;
	int64_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
	
	a0 = A->t[0];
	a1 = A->t[1];
	a2 = A->t[2];
	a3 = A->t[3];
	a4 = A->t[4];
	a5 = A->t[5];
	a6 = A->t[6];
	a7 = A->t[7];
	a8 = A->t[8];
	
	b0 = B->t[1] * LAMBDA;
	b1 = B->t[2] * LAMBDA;
	b2 = B->t[3] * LAMBDA;
	b3 = B->t[4] * LAMBDA;
	b4 = B->t[5] * LAMBDA;
	b5 = B->t[6] * LAMBDA;
	b6 = B->t[7] * LAMBDA;
	b7 = B->t[8] * LAMBDA;
	b8 = B->t[0] * ALPHA;
	b9 = B->t[1] * ALPHA;
	b10 = B->t[2] * ALPHA;
	b11 = B->t[3] * ALPHA;
	b12 = B->t[4] * ALPHA;
	b13 = B->t[5] * ALPHA;
	b14 = B->t[6] * ALPHA;
	b15 = B->t[7] * ALPHA;
	b16 = B->t[8] * ALPHA;
	
	v1m2_0 = a3 - a6;
	v1m2_1 = a4 - a7;
	v1m2_2 = a5 - a8;
	v0m1_0 = a0 - a3;
	v0m1_1 = a1 - a4;
	v0m1_2 = a2 - a5;
	v0m2_0 = a0 - a6;
	v0m2_1 = a1 - a7;
	v0m2_2 = a2 - a8;
	
	m034_0 = b6 + b3 + b0;
	m013_0 = b6 + b3 + b9;
	m012_0 = b6 + b9 + b12;
	m034_1 = b7 + b4 + b1;
	m013_1 = b7 + b4 + b10;
	m012_1 = b7 + b10 + b13;
	m034_2 = b8 + b5 + b2;
	m013_2 = b8 + b5 + b11;
	m012_2 = b8 + b11 + b14;
	m034_3 = b9 + b6 + b3;
	m013_3 = b9 + b6 + b12;
	m012_3 = b9 + b12 + b15;
	m034_4 = b10 + b7 + b4;
	m013_4 = b10 + b7 + b13;
	m012_4 = b10 + b13 + b16;
	
	t3_0 = (__int128) v1m2_0 * b5 + (__int128) v1m2_1 * b4 + (__int128) v1m2_2 * b3;
	t3_1 = (__int128) v1m2_0 * b6 + (__int128) v1m2_1 * b5 + (__int128) v1m2_2 * b4;
	t3_2 = (__int128) v1m2_0 * b7 + (__int128) v1m2_1 * b6 + (__int128) v1m2_2 * b5;
	t4_0 = (__int128) v0m2_0 * b8 + (__int128) v0m2_1 * b7 + (__int128) v0m2_2 * b6;
	t4_1 = (__int128) v0m2_0 * b9 + (__int128) v0m2_1 * b8 + (__int128) v0m2_2 * b7;
	t4_2 = (__int128) v0m2_0 * b10 + (__int128) v0m2_1 * b9 + (__int128) v0m2_2 * b8;
	t5_0 = (__int128) v0m1_0 * b11 + (__int128) v0m1_1 * b10 + (__int128) v0m1_2 * b9;
	t5_1 = (__int128) v0m1_0 * b12 + (__int128) v0m1_1 * b11 + (__int128) v0m1_2 * b10;
	t5_2 = (__int128) v0m1_0 * b13 + (__int128) v0m1_1 * b12 + (__int128) v0m1_2 * b11;
	
	r0 = (__int128) a6 * m034_2 + (__int128) a7 * m034_1 + (__int128) a8 * m034_0 + t3_0 + t4_0;
	r1 = (__int128) a6 * m034_3 + (__int128) a7 * m034_2 + (__int128) a8 * m034_1 + t3_1 + t4_1;
	r2 = (__int128) a6 * m034_4 + (__int128) a7 * m034_3 + (__int128) a8 * m034_2 + t3_2 + t4_2;
	r3 = (__int128) a3 * m013_2 + (__int128) a4 * m013_1 + (__int128) a5 * m013_0 + t5_0 - t3_0;
	r4 = (__int128) a3 * m013_3 + (__int128) a4 * m013_2 + (__int128) a5 * m013_1 + t5_1 - t3_1;
	r5 = (__int128) a3 * m013_4 + (__int128) a4 * m013_3 + (__int128) a5 * m013_2 + t5_2 - t3_2;
	r6 = (__int128) a0 * m012_2 + (__int128) a1 * m012_1 + (__int128) a2 * m012_0 - t4_0 - t5_0;
	r7 = (__int128) a0 * m012_3 + (__int128) a1 * m012_2 + (__int128) a2 * m012_1 - t4_1 - t5_1;
	r8 = (__int128) a0 * m012_4 + (__int128) a1 * m012_3 + (__int128) a2 * m012_2 - t4_2 - t5_2;
	
	t0 = ((uint64_t)r1 + (uint64_t)((uint64_t)r2*GAMMA));
	t1 = ((uint64_t)r2 + (uint64_t)((uint64_t)r3*GAMMA));
	t2 = ((uint64_t)r3 + (uint64_t)((uint64_t)r4*GAMMA));
	t3 = ((uint64_t)r4 + (uint64_t)((uint64_t)r5*GAMMA));
	t4 = ((uint64_t)r5 + (uint64_t)((uint64_t)r6*GAMMA));
	t5 = ((uint64_t)r6 + (uint64_t)((uint64_t)r7*GAMMA));
	t6 = ((uint64_t)r7 + (uint64_t)((uint64_t)r8*GAMMA));
	t7 = ((uint64_t)r8 + (uint64_t)((uint64_t)r0 * GAMMALAMM1));
	t8 = ((uint64_t)((uint64_t)r0 * ONELAMM1) + (uint64_t)((uint64_t)r1 * GAMMALAMM2));
	
	r0 += (__int128) t0 * GAMMA - (__int128)t8 * LAMED;
	r1 += (__int128) t1 * GAMMA - t0;
	r2 += (__int128) t2 * GAMMA - t1;
	r3 += (__int128) t3 * GAMMA - t2;
	r4 += (__int128) t4 * GAMMA - t3;
	r5 += (__int128) t5 * GAMMA - t4;
	r6 += (__int128) t6 * GAMMA - t5;
	r7 += (__int128) t7 * GAMMA - t6;
	r8 += (__int128) t8 * GIMEL - t7;
	
	res->t[0] = (int64_t)(r0>>64);
	res->t[1] = (int64_t)(r1>>64);
	res->t[2] = (int64_t)(r2>>64);
	res->t[3] = (int64_t)(r3>>64);
	res->t[4] = (int64_t)(r4>>64);
	res->t[5] = (int64_t)(r5>>64);
	res->t[6] = (int64_t)(r6>>64);
	res->t[7] = (int64_t)(r7>>64);
	res->t[8] = (int64_t)(r8>>64);
}

void pmns_linearred_mult9(restrict poly res, const restrict poly A, const restrict poly B)
{
	_pmns_linearred_nine(res, A, B);
}

void pmns_linearred_multa9(restrict poly res, const restrict poly A, const restrict poly B)
{
	_pmns_bino_linear_nine(res, A, B);
}

void pmns_linearred_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	switch(N)
	{
		case 4:
			_pmns_linearred_mult(4);
			break;
		case 5:
			_pmns_linearred_mult(5);
			break;
		case 6:
			_pmns_linearred_mult(6)
			break;
		case 7:
			_pmns_linearred_mult(7)
			break;
		case 8:
			_pmns_linearred_toep(8)
			break;
		case 9:
			_pmns_linearred_nine(res, A, B);
			break;
		#ifdef VARIOUSNVALUES
		case 10:
			_pmns_linearred_toep(10)
			break;
		case 11:
			_pmns_linearred_mult(11)
			break;
		case 12:
			_pmns_linearred_toep(12)
			break;
		case 13:
			_pmns_linearred_mult(13)
			break;
		case 14:
			_pmns_linearred_toep(14)
			break;
		#endif
		#ifdef LARGEPMNS
		case 18:
			_pmns_linearred_toep(18)
			break;
		case 36:
			_pmns_linearred_toep(36)
			break;
		case 72:
			_pmns_linearred_toep(72)
			break;
		case 144:
			_pmns_linearred_toep(144)
			break;
		#endif
		default:
			break;
	}
}

void pmns_doublesparse_multa9(restrict poly res, const restrict poly A, const restrict poly B)
{
	_pmns_bino_sparse_nine(res, A, B);
}

void pmns_doublesparse_mult9(restrict poly res, const restrict poly A, const restrict poly B)
{
	_pmns_doublesparse_nine(res, A, B);
}

void pmns_doublesparse_mult(restrict poly res, const restrict poly A, const restrict poly B)
{
	switch(N)
	{
		case 4:
			_pmns_doublesparse_mult(4);
			break;
		case 5:
			_pmns_doublesparse_five(res, A, B);
			break;
		case 6:
			_pmns_doublesparse_mult(6);
			break;
		case 7:
			_pmns_doublesparse_mult(7)
			break;
		case 8:
			_pmns_doublesparse_toep(8);
			break;
		case 9:
			_pmns_doublesparse_nine(res, A, B);
			break;
		#ifdef VARIOUSNVALUES
		case 10:
			_pmns_doublesparse_toep(10)
			break;
		case 11:
			_pmns_doublesparse_mult(11)
			break;
		case 12:
			_pmns_doublesparse_toep(12)
			break;
		case 13:
			_pmns_doublesparse_mult(13)
			break;
		case 14:
			_pmns_doublesparse_toep(14)
			break;
		#endif
		#ifdef LARGEPMNS
		case 18:
			_pmns_doublesparse_toep(18)
			break;
		case 36:
			_pmns_doublesparse_toep(36)
			break;
		case 72:
			_pmns_doublesparse_toep(72)
			break;
		case 144:
			_pmns_doublesparse_toep(144)
			break;
		#endif
		default:
			break;
	}
}

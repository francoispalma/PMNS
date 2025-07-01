#ifndef AVX512TOEPLITZMACROS_H_INCLUDED
#define AVX512TOEPLITZMACROS_H_INCLUDED

void _mm512_madd52_epu64(__m512i* rethi, __m512i* retlo, __m512i op1, __m512i op2)
{
	*retlo = _mm512_madd52lo_epu64(*retlo, op1, op2);
	*rethi = _mm512_madd52hi_epu64(*rethi, op1, op2);
}

void _mm512_madd52_epi64(__m512i* rethi, __m512i* retlo, __m512i op1, __m512i op2)
{
	*retlo = _mm512_madd52lo_epu64(*retlo, op1, op2);
	*rethi = _mm512_madd52hi_epu64(*rethi, op1, op2);
	*rethi = _mm512_mask_sub_epi64(*rethi, _mm512_movepi64_mask(op2), *rethi, op1);
}

void v_mm512_madd52_epi64(__m512i* rethi, __m512i* retlo, __m512i op1, __m512i op2)
{
	*retlo = _mm512_madd52lo_epu64(*retlo, op1, op2);
	*rethi = _mm512_madd52hi_epu64(*rethi, op1, op2);
	*rethi = _mm512_mask_sub_epi64(*rethi, _mm512_movepi64_mask(op1), *rethi, op2);
}

#define TOEP22TOP(X, F, G) {\
	__m512i t0[X/16], t1[X/16], t2[X/16];\
	__m512i t0hi[X/16], t1hi[X/16], t2hi[X/16];\
	__m512i v0p1[X/16], m0m1[X/8], m0m2[X/8];\
	for(int i = 0; i < X/16; i++)\
	{\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[i + X/16]);\
	}\
	for(int i = 0; i < X/8; i++)\
	{\
		m0m1[i] = _mm512_sub_epi64(matr[i+X/16],matr[i+X/8]);\
		m0m2[i] = _mm512_sub_epi64(matr[i+X/16],matr[i]);\
	}\
	F (t0hi, t0, v0p1, matr + X/16);\
	G (t1hi, t1, vect, m0m1);\
	G (t2hi, t2, vect + X/16, m0m2);\
	for(int i = 0; i < X/16; i++)\
	{\
		roplo[i] = _mm512_sub_epi64(t0[i],t2[i]);\
		rophi[i] = _mm512_sub_epi64(t0hi[i],t2hi[i]);\
		roplo[i + X/16] = _mm512_sub_epi64(t0[i],t1[i]);\
		rophi[i + X/16] = _mm512_sub_epi64(t0hi[i],t1hi[i]);\
	}\
}

#define pTOEP22TOP(X, F, G) {\
	__m512i t0[X/16], t1[X/16], t2[X/16];\
	__m512i t0hi[X/16], t1hi[X/16], t2hi[X/16];\
	__m512i v0p1[X/16], m0m1[X/8], m0m2[X/8];\
	for(int i = 0; i < X/16; i++)\
	{\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[i + X/16]);\
	}\
	for(int i = 0; i < X/8; i++)\
	{\
		m0m1[i] = _mm512_sub_epi64(matr[i+X/16],matr[i+X/8]);\
		m0m2[i] = _mm512_sub_epi64(matr[i+X/16],matr[i]);\
	}\
	F (t0hi, t0, v0p1, matr + X/16);\
	G (t1hi, t1, vect, m0m1);\
	G (t2hi, t2, vect + X/16, m0m2);\
	for(int i = 0; i < X/16; i++)\
	{\
		roplo[i] = _mm512_add_epi64(roplo[i],_mm512_sub_epi64(t0[i],t2[i]));\
		rophi[i] = _mm512_add_epi64(rophi[i],_mm512_sub_epi64(t0hi[i],t2hi[i]));\
		roplo[i + X/16] = _mm512_add_epi64(roplo[i + X/16],_mm512_sub_epi64(t0[i],t1[i]));\
		rophi[i + X/16] = _mm512_add_epi64(rophi[i + X/16],_mm512_sub_epi64(t0hi[i],t1hi[i]));\
	}\
}

#define vTOEP22TOP(X, F, G) {\
	__m512i t0[X/16], t1[X/16], t2[X/16];\
	__m512i t0hi[X/16], t1hi[X/16], t2hi[X/16];\
	__m512i v0p1[X/16], m0m1[X/8], m0m2[X/8];\
	for(int i = 0; i < X/16; i++)\
	{\
		v0p1[i] = _mm512_sub_epi64(vect[i],vect[i + X/16]);\
	}\
	for(int i = 0; i < X/8; i++)\
	{\
		m0m1[i] = _mm512_add_epi64(matr[i+X/16],matr[i+X/8]);\
		m0m2[i] = _mm512_add_epi64(matr[i+X/16],matr[i]);\
	}\
	F (t0hi, t0, v0p1, matr + X/16);\
	G (t1hi, t1, vect, m0m1);\
	G (t2hi, t2, vect + X/16, m0m2);\
	for(int i = 0; i < X/16; i++)\
	{\
		roplo[i] = _mm512_add_epi64(t0[i],t2[i]);\
		rophi[i] = _mm512_add_epi64(t0hi[i],t2hi[i]);\
		roplo[i + X/16] = _mm512_sub_epi64(t1[i],t0[i]);\
		rophi[i + X/16] = _mm512_sub_epi64(t1hi[i],t0hi[i]);\
	}\
}

#define vpTOEP22TOP(X, F, G) {\
	__m512i t0[X/16], t1[X/16], t2[X/16];\
	__m512i t0hi[X/16], t1hi[X/16], t2hi[X/16];\
	__m512i v0p1[X/16], m0m1[X/8], m0m2[X/8];\
	for(int i = 0; i < X/16; i++)\
	{\
		v0p1[i] = _mm512_sub_epi64(vect[i],vect[i + X/16]);\
	}\
	for(int i = 0; i < X/8; i++)\
	{\
		m0m1[i] = _mm512_add_epi64(matr[i+X/16],matr[i+X/8]);\
		m0m2[i] = _mm512_add_epi64(matr[i+X/16],matr[i]);\
	}\
	F (t0hi, t0, v0p1, matr + X/16);\
	G (t1hi, t1, vect, m0m1);\
	G (t2hi, t2, vect + X/16, m0m2);\
	for(int i = 0; i < X/16; i++)\
	{\
		roplo[i] = _mm512_add_epi64(roplo[i],_mm512_add_epi64(t0[i],t2[i]));\
		rophi[i] = _mm512_add_epi64(rophi[i],_mm512_add_epi64(t0hi[i],t2hi[i]));\
		roplo[i + X/16] = _mm512_add_epi64(roplo[i + X/16],_mm512_sub_epi64(t1[i],t0[i]));\
		rophi[i + X/16] = _mm512_add_epi64(rophi[i + X/16],_mm512_sub_epi64(t1hi[i],t0hi[i]));\
	}\
}

#define TOEP33TOP(X, F, G) {\
	__m512i t0[X/24], t1[X/24], t2[X/24], t3[X/24], t4[X/24], t5[X/24];\
	__m512i t0hi[X/24], t1hi[X/24], t2hi[X/24], t3hi[X/24], t4hi[X/24], t5hi[X/24];\
	__m512i v1p2[X/24], v0p2[X/24], v0p1[X/24], m034[X/12], m013[X/12], m012[X/12];\
	for(int i = 0; i < X/24; i++)\
	{\
		v1p2[i] = _mm512_add_epi64(vect[X/24 + i],vect[2*X/24 + i]);\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[X/24 + i]);\
		v0p2[i] = _mm512_add_epi64(vect[i],vect[2*X/24 + i]);\
	}\
	for(int i = 0; i < X/12; i++)\
	{\
		m034[i] = _mm512_sub_epi64(matr[i],_mm512_add_epi64(matr[i + 2*X/24],matr[i + X/24]));\
		m013[i] = _mm512_sub_epi64(matr[i+2*X/24],_mm512_add_epi64(matr[i + X/8],matr[i + X/24]));\
		m012[i] = _mm512_sub_epi64(matr[i + 4*X/24],_mm512_add_epi64(matr[i + 2*X/24],matr[i + X/8]));\
	}\
	G (t0hi, t0, vect + 2*X/24, m034);\
	G (t1hi, t1, vect + X/24, m013);\
	G (t2hi, t2, vect, m012);\
	F (t3hi, t3, v1p2,  matr + X/24);\
	F (t4hi, t4, v0p2, matr + 2*X/24);\
	F (t5hi, t5, v0p1,  matr + X/8);\
	for(int i = 0; i < X/24; i++)\
	{\
		roplo[i] = _mm512_add_epi64(t0[i],_mm512_add_epi64(t3[i],t4[i]));\
		rophi[i] = _mm512_add_epi64(t0hi[i],_mm512_add_epi64(t3hi[i],t4hi[i]));\
		roplo[i+X/24] = _mm512_add_epi64(t1[i],_mm512_add_epi64(t3[i],t5[i]));\
		rophi[i+X/24] = _mm512_add_epi64(t1hi[i],_mm512_add_epi64(t3hi[i],t5hi[i]));\
		roplo[i+2*X/24] = _mm512_add_epi64(t2[i],_mm512_add_epi64(t4[i],t5[i]));\
		rophi[i+2*X/24] = _mm512_add_epi64(t2hi[i],_mm512_add_epi64(t4hi[i],t5hi[i]));\
	}\
}

#define pTOEP33TOP(X, F, G) {\
	__m512i t0[X/24], t1[X/24], t2[X/24], t3[X/24], t4[X/24], t5[X/24];\
	__m512i t0hi[X/24], t1hi[X/24], t2hi[X/24], t3hi[X/24], t4hi[X/24], t5hi[X/24];\
	__m512i v1p2[X/24], v0p2[X/24], v0p1[X/24], m034[X/12], m013[X/12], m012[X/12];\
	for(int i = 0; i < X/24; i++)\
	{\
		v1p2[i] = _mm512_add_epi64(vect[X/24 + i],vect[2*X/24 + i]);\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[X/24 + i]);\
		v0p2[i] = _mm512_add_epi64(vect[i],vect[2*X/24 + i]);\
	}\
	for(int i = 0; i < X/12; i++)\
	{\
		m034[i] = _mm512_sub_epi64(matr[i],_mm512_add_epi64(matr[i + 2*X/24],matr[i + X/24]));\
		m013[i] = _mm512_sub_epi64(matr[i+2*X/24],_mm512_add_epi64(matr[i + X/8],matr[i + X/24]));\
		m012[i] = _mm512_sub_epi64(matr[i + 4*X/24],_mm512_add_epi64(matr[i + 2*X/24],matr[i + X/8]));\
	}\
	G (t0hi, t0, vect + 2*X/24, m034);\
	G (t1hi, t1, vect + X/24, m013);\
	G (t2hi, t2, vect, m012);\
	F (t3hi, t3, v1p2,  matr + X/24);\
	F (t4hi, t4, v0p2, matr + 2*X/24);\
	F (t5hi, t5, v0p1,  matr + X/8);\
	for(int i = 0; i < X/24; i++)\
	{\
		roplo[i] = _mm512_add_epi64(roplo[i],_mm512_add_epi64(t0[i],_mm512_add_epi64(t3[i],t4[i])));\
		rophi[i] = _mm512_add_epi64(rophi[i],_mm512_add_epi64(t0hi[i],_mm512_add_epi64(t3hi[i],t4hi[i])));\
		roplo[i+X/24] = _mm512_add_epi64(roplo[i+X/24],_mm512_add_epi64(t1[i],_mm512_add_epi64(t3[i],t5[i])));\
		rophi[i+X/24] = _mm512_add_epi64(rophi[i+X/24],_mm512_add_epi64(t1hi[i],_mm512_add_epi64(t3hi[i],t5hi[i])));\
		roplo[i+2*X/24] = _mm512_add_epi64(roplo[i+2*X/24],_mm512_add_epi64(t2[i],_mm512_add_epi64(t4[i],t5[i])));\
		rophi[i+2*X/24] = _mm512_add_epi64(rophi[i+2*X/24],_mm512_add_epi64(t2hi[i],_mm512_add_epi64(t4hi[i],t5hi[i])));\
	}\
}

#define vTOEP33TOP(X, F, G) {\
	__m512i t0[X/24], t1[X/24], t2[X/24], t3[X/24], t4[X/24], t5[X/24];\
	__m512i t0hi[X/24], t1hi[X/24], t2hi[X/24], t3hi[X/24], t4hi[X/24], t5hi[X/24];\
	__m512i v1m2[X/24], v0m2[X/24], v0m1[X/24], m034[X/12], m013[X/12], m012[X/12], m03;\
	for(int i = 0; i < X/24; i++)\
	{\
		v1m2[i] = _mm512_sub_epi64(vect[X/24 + i],vect[2*X/24 + i]);\
		v0m1[i] = _mm512_sub_epi64(vect[i],vect[X/24 + i]);\
		v0m2[i] = _mm512_sub_epi64(vect[i],vect[2*X/24 + i]);\
	}\
	for(int i = 0; i < X/12; i++)\
	{\
		m03 = _mm512_add_epi64(matr[i + 2*X/24],matr[i + X/24]);\
		m034[i] = _mm512_add_epi64(matr[i], m03);\
		m013[i] = _mm512_add_epi64(matr[i + 3*X/24], m03);\
		m012[i] = _mm512_add_epi64(matr[i + 4*X/24], _mm512_add_epi64(matr[i + 2*X/24],matr[i + 3*X/24]));\
	}\
	G (t0hi, t0, vect + 2*X/24, m034);\
	G (t1hi, t1, vect + X/24, m013);\
	G (t2hi, t2, vect, m012);\
	F (t3hi, t3, v1m2,  matr + X/24);\
	F (t4hi, t4, v0m2, matr + 2*X/24);\
	F (t5hi, t5, v0m1,  matr + X/8);\
	for(int i = 0; i < X/24; i++)\
	{\
		roplo[i] = _mm512_add_epi64(t0[i],_mm512_add_epi64(t3[i],t4[i]));\
		rophi[i] = _mm512_add_epi64(t0hi[i],_mm512_add_epi64(t3hi[i],t4hi[i]));\
		roplo[i+X/24] = _mm512_add_epi64(t1[i],_mm512_sub_epi64(t5[i],t3[i]));\
		rophi[i+X/24] = _mm512_add_epi64(t1hi[i],_mm512_sub_epi64(t5hi[i],t3hi[i]));\
		roplo[i+2*X/24] = _mm512_sub_epi64(t2[i],_mm512_add_epi64(t4[i],t5[i]));\
		rophi[i+2*X/24] = _mm512_sub_epi64(t2hi[i],_mm512_add_epi64(t4hi[i],t5hi[i]));\
	}\
}

#define vpTOEP33TOP(X, F, G) {\
	__m512i t0[X/24], t1[X/24], t2[X/24], t3[X/24], t4[X/24], t5[X/24];\
	__m512i t0hi[X/24], t1hi[X/24], t2hi[X/24], t3hi[X/24], t4hi[X/24], t5hi[X/24];\
	__m512i v1m2[X/24], v0m2[X/24], v0m1[X/24], m034[X/12], m013[X/12], m012[X/12], m03;\
	for(int i = 0; i < X/24; i++)\
	{\
		v1m2[i] = _mm512_sub_epi64(vect[X/24 + i],vect[2*X/24 + i]);\
		v0m1[i] = _mm512_sub_epi64(vect[i],vect[X/24 + i]);\
		v0m2[i] = _mm512_sub_epi64(vect[i],vect[2*X/24 + i]);\
	}\
	for(int i = 0; i < X/12; i++)\
	{\
		m03 = _mm512_add_epi64(matr[i + 2*X/24],matr[i + X/24]);\
		m034[i] = _mm512_add_epi64(matr[i], m03);\
		m013[i] = _mm512_add_epi64(matr[i + X/8], m03);\
		m012[i] = _mm512_add_epi64(matr[i + 4*X/24], _mm512_add_epi64(matr[i + 2*X/24],matr[i + X/8]));\
	}\
	G (t0hi, t0, vect + 2*X/24, m034);\
	G (t1hi, t1, vect + X/24, m013);\
	G (t2hi, t2, vect, m012);\
	F (t3hi, t3, v1m2,  matr + X/24);\
	F (t4hi, t4, v0m2, matr + 2*X/24);\
	F (t5hi, t5, v0m1,  matr + X/8);\
	for(int i = 0; i < X/24; i++)\
	{\
		roplo[i] = _mm512_add_epi64(roplo[i],_mm512_add_epi64(t0[i],_mm512_add_epi64(t3[i],t4[i])));\
		rophi[i] = _mm512_add_epi64(rophi[i],_mm512_add_epi64(t0hi[i],_mm512_add_epi64(t3hi[i],t4hi[i])));\
		roplo[i+X/24] = _mm512_add_epi64(roplo[i+X/24],_mm512_add_epi64(t1[i],_mm512_sub_epi64(t5[i],t3[i])));\
		rophi[i+X/24] = _mm512_add_epi64(rophi[i+X/24],_mm512_add_epi64(t1hi[i],_mm512_sub_epi64(t5hi[i],t3hi[i])));\
		roplo[i+2*X/24] = _mm512_add_epi64(roplo[i+2*X/24],_mm512_sub_epi64(t2[i],_mm512_add_epi64(t4[i],t5[i])));\
		rophi[i+2*X/24] = _mm512_add_epi64(rophi[i+2*X/24],_mm512_sub_epi64(t2hi[i],_mm512_add_epi64(t4hi[i],t5hi[i])));\
	}\
}

#define TOEP55TOP(X, F, G) {\
	__m512i t0[X/40], t1[X/40], t2[X/40], t3[X/40], t4[X/40], t5[X/40], t6[X/40], t7[X/40], t8[X/40], t9[X/40], t10[X/40], t11[X/40], t12[X/40], t13[X/40], t14[X/40];\
	__m512i t0hi[X/40], t1hi[X/40], t2hi[X/40], t3hi[X/40], t4hi[X/40], t5hi[X/40], t6hi[X/40], t7hi[X/40], t8hi[X/40], t9hi[X/40], t10hi[X/40], t11hi[X/40], t12hi[X/40], t13hi[X/40], t14hi[X/40];\
	__m512i v0p4[X/40], v1p4[X/40], v2p4[X/40], v3p4[X/40], v0p3[X/40], v1p3[X/40], v2p3[X/40], v0p2[X/40], v1p2[X/40], v0p1[X/40];\
	__m512i m01234[X/20], m50123[X/20], m65012[X/20], m76501[X/20], m87650[X/20];\
	__m512i m50, m750, m12;\
	for(int i = 0; i < X/40; i++)\
	{\
		v0p4[i] = _mm512_add_epi64(vect[i],vect[i+4*X/40]);\
		v1p4[i] = _mm512_add_epi64(vect[i+X/40],vect[i+4*X/40]);\
		v2p4[i] = _mm512_add_epi64(vect[i+2*X/40],vect[i+4*X/40]);\
		v3p4[i] = _mm512_add_epi64(vect[i+3*X/40],vect[i+4*X/40]);\
		v0p3[i] = _mm512_add_epi64(vect[i],vect[i+3*X/40]);\
		v1p3[i] = _mm512_add_epi64(vect[i+X/40],vect[i+3*X/40]);\
		v2p3[i] = _mm512_add_epi64(vect[i+2*X/40],vect[i+3*X/40]);\
		v0p2[i] = _mm512_add_epi64(vect[i],vect[i+2*X/40]);\
		v1p2[i] = _mm512_add_epi64(vect[i+X/40],vect[i+2*X/40]);\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[i+X/40]);\
	}\
	for(int i = 0; i < X/20; i++)\
	{\
		m50 = _mm512_add_epi64(matr[i+3*X/40],matr[i+4*X/40]);\
		m750 = _mm512_add_epi64(matr[i+1*X/40],m50);\
		m12 = _mm512_add_epi64(matr[i+5*X/40],matr[i+6*X/40]);\
		m01234[i] = _mm512_sub_epi64(matr[i+8*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+4*X/40],m12),matr[i+7*X/40]));\
		m50123[i] = _mm512_sub_epi64(matr[i+6*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+5*X/40],m50),matr[i+7*X/40]));\
		m65012[i] = _mm512_sub_epi64(matr[i+4*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+3*X/40],m12),matr[i+2*X/40]));\
		m76501[i] = _mm512_sub_epi64(matr[i+2*X/40],_mm512_add_epi64(matr[i+5*X/40],m750));\
		m87650[i] = _mm512_sub_epi64(matr[i],_mm512_add_epi64(matr[i+2*X/40],m750));\
	}\
	G (t0hi, t0, vect + 4*X/40, m87650);\
	G (t1hi, t1, vect + 3*X/40, m76501);\
	G (t2hi, t2, vect + 2*X/40, m65012);\
	G (t3hi, t3, vect + 1*X/40, m50123);\
	G (t4hi, t4, vect        , m01234);\
	F (t5hi, t5 , v0p4, matr + 4*X/40);\
	F (t6hi, t6 , v1p4, matr + 3*X/40);\
	F (t7hi, t7 , v2p4, matr + 2*X/40);\
	F (t8hi, t8 , v3p4, matr + 1*X/40);\
	F (t9hi, t9 , v0p3, matr + 5*X/40);\
	F (t10hi, t10, v1p3, matr + 4*X/40);\
	F (t11hi, t11, v2p3, matr + 3*X/40);\
	F (t12hi, t12, v0p2, matr + 6*X/40);\
	F (t13hi, t13, v1p2, matr + 5*X/40);\
	F (t14hi, t14, v0p1, matr + 7*X/40);\
	for(int i = 0; i < X/40; i++)\
	{\
		roplo[i] = _mm512_add_epi64(t0[i], _mm512_add_epi64(t5[i],_mm512_add_epi64(t6[i],_mm512_add_epi64(t7[i],t8[i]))));\
		rophi[i] = _mm512_add_epi64(t0hi[i], _mm512_add_epi64(t5hi[i],_mm512_add_epi64(t6hi[i],_mm512_add_epi64(t7hi[i],t8hi[i]))));\
		roplo[i+X/40] = _mm512_add_epi64(t1[i], _mm512_add_epi64(t8[i],_mm512_add_epi64(t9[i],_mm512_add_epi64(t10[i],t11[i]))));\
		rophi[i+X/40] = _mm512_add_epi64(t1hi[i], _mm512_add_epi64(t8hi[i],_mm512_add_epi64(t9hi[i],_mm512_add_epi64(t10hi[i],t11hi[i]))));\
		roplo[i+2*X/40] = _mm512_add_epi64(t2[i], _mm512_add_epi64(t7[i],_mm512_add_epi64(t11[i],_mm512_add_epi64(t12[i],t13[i]))));\
		rophi[i+2*X/40] = _mm512_add_epi64(t2hi[i], _mm512_add_epi64(t7hi[i],_mm512_add_epi64(t11hi[i],_mm512_add_epi64(t12hi[i],t13hi[i]))));\
		roplo[i+3*X/40] = _mm512_add_epi64(t3[i], _mm512_add_epi64(t6[i],_mm512_add_epi64(t10[i],_mm512_add_epi64(t13[i],t14[i]))));\
		rophi[i+3*X/40] = _mm512_add_epi64(t3hi[i], _mm512_add_epi64(t6hi[i],_mm512_add_epi64(t10hi[i],_mm512_add_epi64(t13hi[i],t14hi[i]))));\
		roplo[i+4*X/40] = _mm512_add_epi64(t4[i], _mm512_add_epi64(t5[i],_mm512_add_epi64(t9[i],_mm512_add_epi64(t12[i],t14[i]))));\
		rophi[i+4*X/40] = _mm512_add_epi64(t4hi[i], _mm512_add_epi64(t5hi[i],_mm512_add_epi64(t9hi[i],_mm512_add_epi64(t12hi[i],t14hi[i]))));\
	}\
}

#define pTOEP55TOP(X, F, G) {\
	__m512i t0[X/8], t1[X/8], t2[X/8], t3[X/8], t4[X/8], t5[X/8], t6[X/8], t7[X/8], t8[X/8], t9[X/8], t10[X/8], t11[X/8], t12[X/8], t13[X/8], t14[X/8];\
	__m512i t0hi[X/8], t1hi[X/8], t2hi[X/8], t3hi[X/8], t4hi[X/8], t5hi[X/8], t6hi[X/8], t7hi[X/8], t8hi[X/8], t9hi[X/8], t10hi[X/8], t11hi[X/8], t12hi[X/8], t13hi[X/8], t14hi[X/8];\
	__m512i v0p4[X/40], v1p4[X/40], v2p4[X/40], v3p4[X/40], v0p3[X/40], v1p3[X/40], v2p3[X/40], v0p2[X/40], v1p2[X/40], v0p1[X/40];\
	__m512i m01234[X/20], m50123[X/20], m65012[X/20], m76501[X/20], m87650[X/20];\
	__m512i m50, m750, m12;\
	for(int i = 0; i < X/40; i++)\
	{\
		v0p4[i] = _mm512_add_epi64(vect[i],vect[i+4*X/40]);\
		v1p4[i] = _mm512_add_epi64(vect[i+X/40],vect[i+4*X/40]);\
		v2p4[i] = _mm512_add_epi64(vect[i+2*X/40],vect[i+4*X/40]);\
		v3p4[i] = _mm512_add_epi64(vect[i+3*X/40],vect[i+4*X/40]);\
		v0p3[i] = _mm512_add_epi64(vect[i],vect[i+3*X/40]);\
		v1p3[i] = _mm512_add_epi64(vect[i+X/40],vect[i+3*X/40]);\
		v2p3[i] = _mm512_add_epi64(vect[i+2*X/40],vect[i+3*X/40]);\
		v0p2[i] = _mm512_add_epi64(vect[i],vect[i+2*X/40]);\
		v1p2[i] = _mm512_add_epi64(vect[i+X/40],vect[i+2*X/40]);\
		v0p1[i] = _mm512_add_epi64(vect[i],vect[i+X/40]);\
	}\
	for(int i = 0; i < X/20; i++)\
	{\
		m50 = _mm512_add_epi64(matr[i+3*X/40],matr[i+4*X/40]);\
		m750 = _mm512_add_epi64(matr[i+1*X/40],m50);\
		m12 = _mm512_add_epi64(matr[i+5*X/40],matr[i+6*X/40]);\
		m01234[i] = _mm512_sub_epi64(matr[i+8*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+4*X/40],m12),matr[i+7*X/40]));\
		m50123[i] = _mm512_sub_epi64(matr[i+6*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+5*X/40],m50),matr[i+7*X/40]));\
		m65012[i] = _mm512_sub_epi64(matr[i+4*X/40],_mm512_add_epi64(_mm512_add_epi64(matr[i+3*X/40],m12),matr[i+2*X/40]));\
		m76501[i] = _mm512_sub_epi64(matr[i+2*X/40],_mm512_add_epi64(matr[i+5*X/40],m750));\
		m87650[i] = _mm512_sub_epi64(matr[i],_mm512_add_epi64(matr[i+2*X/40],m750));\
	}\
	G (t0hi, t0, vect + 4*X/40, m87650);\
	G (t1hi, t1, vect + 3*X/40, m76501);\
	G (t2hi, t2, vect + 2*X/40, m65012);\
	G (t3hi, t3, vect + 1*X/40, m50123);\
	G (t4hi, t4, vect        , m01234);\
	F (t5hi, t5 , v0p4, matr + 4*X/40);\
	F (t6hi, t6 , v1p4, matr + 3*X/40);\
	F (t7hi, t7 , v2p4, matr + 2*X/40);\
	F (t8hi, t8 , v3p4, matr + 1*X/40);\
	F (t9hi, t9 , v0p3, matr + 5*X/40);\
	F (t10hi, t10, v1p3, matr + 4*X/40);\
	F (t11hi, t11, v2p3, matr + 3*X/40);\
	F (t12hi, t12, v0p2, matr + 6*X/40);\
	F (t13hi, t13, v1p2, matr + 5*X/40);\
	F (t14hi, t14, v0p1, matr + 7*X/40);\
	for(int i = 0; i < X/40; i++)\
	{\
		roplo[i] = _mm512_add_epi64(roplo[i],_mm512_add_epi64(t0[i], _mm512_add_epi64(t5[i],_mm512_add_epi64(t6[i],_mm512_add_epi64(t7[i],t8[i])))));\
		rophi[i] = _mm512_add_epi64(rophi[i],_mm512_add_epi64(t0hi[i], _mm512_add_epi64(t5hi[i],_mm512_add_epi64(t6hi[i],_mm512_add_epi64(t7hi[i],t8hi[i])))));\
		roplo[i+X/40] = _mm512_add_epi64(roplo[i+X/40],_mm512_add_epi64(t1[i], _mm512_add_epi64(t8[i],_mm512_add_epi64(t9[i],_mm512_add_epi64(t10[i],t11[i])))));\
		rophi[i+X/40] = _mm512_add_epi64(rophi[i+X/40],_mm512_add_epi64(t1hi[i], _mm512_add_epi64(t8hi[i],_mm512_add_epi64(t9hi[i],_mm512_add_epi64(t10hi[i],t11hi[i])))));\
		roplo[i+2*X/40] = _mm512_add_epi64(roplo[i+2*X/40],_mm512_add_epi64(t2[i], _mm512_add_epi64(t7[i],_mm512_add_epi64(t11[i],_mm512_add_epi64(t12[i],t13[i])))));\
		rophi[i+2*X/40] = _mm512_add_epi64(rophi[i+2*X/40],_mm512_add_epi64(t2hi[i], _mm512_add_epi64(t7hi[i],_mm512_add_epi64(t11hi[i],_mm512_add_epi64(t12hi[i],t13hi[i])))));\
		roplo[i+3*X/40] = _mm512_add_epi64(roplo[i+3*X/40],_mm512_add_epi64(t3[i], _mm512_add_epi64(t6[i],_mm512_add_epi64(t10[i],_mm512_add_epi64(t13[i],t14[i])))));\
		rophi[i+3*X/40] = _mm512_add_epi64(rophi[i+3*X/40],_mm512_add_epi64(t3hi[i], _mm512_add_epi64(t6hi[i],_mm512_add_epi64(t10hi[i],_mm512_add_epi64(t13hi[i],t14hi[i])))));\
		roplo[i+4*X/40] = _mm512_add_epi64(roplo[i+4*X/40],_mm512_add_epi64(t4[i], _mm512_add_epi64(t5[i],_mm512_add_epi64(t9[i],_mm512_add_epi64(t12[i],t14[i])))));\
		rophi[i+4*X/40] = _mm512_add_epi64(rophi[i+4*X/40],_mm512_add_epi64(t4hi[i], _mm512_add_epi64(t5hi[i],_mm512_add_epi64(t9hi[i],_mm512_add_epi64(t12hi[i],t14hi[i])))));\
	}\
}

#define SCHOOLBOOK(X) {\
	int64_t vecv[X/8][8];\
	for(int i = 0; i < X/8; i++)\
		_mm512_store_si512((__m512i*)vecv[i], vect[i]);\
	for(int i = 0; i < X/8; i++)\
	{\
		roplo[i] = _mm512_setzero_si512();\
		rophi[i] = _mm512_setzero_si512();\
		for(int j = 0; j < X/8; j++)\
		{\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][0]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-0));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][1]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-1));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][2]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-2));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][3]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-3));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][4]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-4));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][5]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-5));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][6]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-6));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][7]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-7));\
		}\
	}\
}\

#define sSCHOOLBOOK(X) {\
	int64_t vecv[X/8][8];\
	for(int i = 0; i < X/8; i++)\
		_mm512_store_si512((__m512i*)vecv[i], vect[i]);\
	for(int i = 0; i < X/8; i++)\
	{\
		roplo[i] = _mm512_setzero_si512();\
		rophi[i] = _mm512_setzero_si512();\
		for(int j = 0; j < X/8; j++)\
		{\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][0]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-0));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][1]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-1));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][2]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-2));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][3]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-3));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][4]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-4));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][5]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-5));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][6]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-6));\
			_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][7]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-7));\
		}\
	}\
}\

#define vSCHOOLBOOK(X) {\
	int64_t vecv[X/8][8];\
	for(int i = 0; i < X/8; i++)\
		_mm512_store_si512((__m512i*)vecv[i], vect[i]);\
	for(int i = 0; i < X/8; i++)\
	{\
		roplo[i] = _mm512_setzero_si512();\
		rophi[i] = _mm512_setzero_si512();\
		for(int j = 0; j < X/8; j++)\
		{\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][0]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-0));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][1]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-1));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][2]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-2));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][3]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-3));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][4]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-4));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][5]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-5));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][6]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-6));\
			v_mm512_madd52_epi64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][7]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-7));\
		}\
	}\
}\

#define pSCHOOLBOOK(X) {\
	int64_t vecv[X/8][8];\
	for(int i = 0; i < X/8; i++)\
		_mm512_store_si512((__m512i*)vecv[i], vect[i]);\
	for(int i = 0; i < X/8; i++)\
	{\
		for(int j = 0; j < X/8; j++)\
		{\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][0]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-0));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][1]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-1));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][2]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-2));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][3]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-3));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][4]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-4));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][5]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-5));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][6]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-6));\
			_mm512_madd52_epu64(rophi + i, roplo + i, _mm512_set1_epi64(vecv[j][7]), _mm512_alignr_epi64(matr[X/8+i-j], matr[X/8+i-1-j],7-7));\
		}\
	}\
}\

#define M1SCHOOLBOOK(X) {\
	__m512i tmp;\
	for(int i = 0; i < X/8; i++)\
		rop[i] = _mm512_setzero_si512();\
	for(int i = 0; i < X; i++)\
	{\
		tmp = _mm512_set1_epi64(((int64_t*)vect)[i]);\
		for(int j = 0; j < X/8; j++)\
		{\
			rop[j] = _mm512_madd52lo_epu64(rop[j], tmp, matr[X - 1 - i + 8*j]);\
		}\
	}\
}

#endif

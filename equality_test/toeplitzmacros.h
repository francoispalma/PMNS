#ifndef TOEPLITZMACROS_H_INCLUDED
#define TOEPLITZMACROS_H_INCLUDED

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

#define TOEP55TOP(X, F) {\
	__int128 t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	int64_t v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	int64_t m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	int64_t m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = vect[i] - vect[4*X/5 + i];\
		v1m4[i] = vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = vect[i] - vect[3*X/5 + i];\
		v1m3[i] = vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = vect[i] - vect[2*X/5 + i];\
		v1m2[i] = vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = m01 + matr[6*X/5 + i];\
		m0123 = m012 + matr[7*X/5 + i];\
		m01234[i] = m0123 + matr[8*X/5 + i];\
		m50123[i] = matr[3*X/5 + i] + m0123;\
		m65 = matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = matr[1*X/5 + i] + m65;\
		m65012[i] = m65 + m012;\
		m76501[i] = m765 + m01;\
		m87650[i] = matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	F (t5 , v0m4, matr + 4*X/5);\
	F (t6 , v1m4, matr + 3*X/5);\
	F (t7 , v2m4, matr + 2*X/5);\
	F (t8 , v3m4, matr + 1*X/5);\
	F (t9 , v0m3, matr + 5*X/5);\
	F (t10, v1m3, matr + 4*X/5);\
	F (t11, v2m3, matr + 3*X/5);\
	F (t12, v0m2, matr + 6*X/5);\
	F (t13, v1m2, matr + 5*X/5);\
	F (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#define MTOEP22TOP(X, F, G) {\
	__int128 t0[X/2], t1[X/2], t2[X/2], v0p1[X/2];\
	int64_t m0m1[X-1], m0m2[X-1];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = (__int128)vect[i] + vect[i + X/2];\
	}\
	for(int i = 0; i < X-1; i++)\
	{\
		m0m1[i] = matr[i + X/2] - matr[i + X];\
		m0m2[i] = matr[i + X/2] - matr[i];\
	}\
	G (t0, v0p1, matr + X/2);\
	F (t1, vect, m0m1);\
	F (t2, vect + X/2, m0m2);\
	for(int i = 0; i < X/2; i++)\
	{\
		rop[i] = t0[i] - t2[i];\
		rop[i + X/2] = t0[i] - t1[i];\
	}\
}

#define MMTOEP22TOP(X, F) {\
	__int128 t0[X/2], t1[X/2], t2[X/2], v0p1[X/2];\
	int64_t m0m1[X-1], m0m2[X-1];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = (__int128)vect[i] + vect[i + X/2];\
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

#define bigmatTOEP22TOP(X, F) {\
	__int128 t0[X/2], t1[X/2], t2[X/2], m0m1[X-1], m0m2[X-1];\
	int64_t v0p1[X/2];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = vect[i] + vect[i + X/2];\
	}\
	for(int i = 0; i < X-1; i++)\
	{\
		m0m1[i] = (__int128)matr[i + X/2] - matr[i + X];\
		m0m2[i] = (__int128)matr[i + X/2] - matr[i];\
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

#define TOEP22TOP128(X, F) {\
	__int128 t0[X/2], t1[X/2], t2[X/2], v0p1[X/2], m0m1[X-1], m0m2[X-1];\
	for(int i = 0; i < X/2; i++)\
	{\
		v0p1[i] = (__int128)vect[i] + vect[i + X/2];\
	}\
	for(int i = 0; i < X-1; i++)\
	{\
		m0m1[i] = (__int128)matr[i + X/2] - matr[i + X];\
		m0m2[i] = (__int128)matr[i + X/2] - matr[i];\
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

#define MTOEP33TOP(X, F, G) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3], v1m2[X/3], v0m2[X/3], v0m1[X/3];\
	int64_t m03, m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = (__int128)vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/3 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/3 + i];\
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
	G (t3, v1m2,  matr + X/3);\
	G (t4, v0m2, matr + 2*X/3);\
	G (t5, v0m1,  matr + X);\
	for(int i = 0; i < X/3; i++)\
	{\
		rop[i] = t0[i] + t3[i] + t4[i];\
		rop[i + X/3] = t1[i] - t3[i] + t5[i];\
		rop[i + 2*X/3] = t2[i] - t4[i] - t5[i];\
	}\
}

#define MMTOEP33TOP(X, F) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3], v1m2[X/3], v0m2[X/3], v0m1[X/3];\
	int64_t m03, m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = (__int128)vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/3 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/3 + i];\
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

#define bigmatTOEP33TOP(X, F) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3], m03, m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1];\
	int64_t v1m2[X/3], v0m2[X/3], v0m1[X/3];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = vect[i] - vect[X/3 + i];\
		v0m2[i] = vect[i] - vect[2*X/3 + i];\
	}\
	for(int i = 0; i < (2*X/3) - 1; i++)\
	{\
		m03 = (__int128)matr[i + 2*X/3] + matr[i + X/3];\
		m034[i] = (__int128)m03 + matr[i];\
		m013[i] = (__int128)m03 +  matr[i + X];\
		m012[i] = (__int128)matr[i + 2*X/3] +  matr[i + X] + matr[i + 4*X/3];\
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

#define TOEP33TOP128(X, F) {\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3], v1m2[X/3], v0m2[X/3], v0m1[X/3], m03, m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1];\
	for(int i = 0; i < X/3; i++)\
	{\
		v1m2[i] = (__int128)vect[X/3 + i] - vect[2*X/3 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/3 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/3 + i];\
	}\
	for(int i = 0; i < (2*X/3) - 1; i++)\
	{\
		m03 = (__int128)matr[i + 2*X/3] + matr[i + X/3];\
		m034[i] = (__int128)m03 + matr[i];\
		m013[i] = (__int128)m03 +  matr[i + X];\
		m012[i] = (__int128)matr[i + 2*X/3] +  matr[i + X] + matr[i + 4*X/3];\
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

#define MTOEP55TOP(X, F, G) {\
	__int128 t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	__int128 v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	int64_t m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	int64_t m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = (__int128)vect[i] - vect[4*X/5 + i];\
		v1m4[i] = (__int128)vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = (__int128)vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = (__int128)vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = (__int128)vect[i] - vect[3*X/5 + i];\
		v1m3[i] = (__int128)vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = (__int128)vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/5 + i];\
		v1m2[i] = (__int128)vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = m01 + matr[6*X/5 + i];\
		m0123 = m012 + matr[7*X/5 + i];\
		m01234[i] = m0123 + matr[8*X/5 + i];\
		m50123[i] = matr[3*X/5 + i] + m0123;\
		m65 = matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = matr[1*X/5 + i] + m65;\
		m65012[i] = m65 + m012;\
		m76501[i] = m765 + m01;\
		m87650[i] = matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	G (t5 , v0m4, matr + 4*X/5);\
	G (t6 , v1m4, matr + 3*X/5);\
	G (t7 , v2m4, matr + 2*X/5);\
	G (t8 , v3m4, matr + 1*X/5);\
	G (t9 , v0m3, matr + 5*X/5);\
	G (t10, v1m3, matr + 4*X/5);\
	G (t11, v2m3, matr + 3*X/5);\
	G (t12, v0m2, matr + 6*X/5);\
	G (t13, v1m2, matr + 5*X/5);\
	G (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#define MMTOEP55TOP(X, F) {\
	__int128 t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	__int128 v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	int64_t m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	int64_t m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = (__int128)vect[i] - vect[4*X/5 + i];\
		v1m4[i] = (__int128)vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = (__int128)vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = (__int128)vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = (__int128)vect[i] - vect[3*X/5 + i];\
		v1m3[i] = (__int128)vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = (__int128)vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/5 + i];\
		v1m2[i] = (__int128)vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = m01 + matr[6*X/5 + i];\
		m0123 = m012 + matr[7*X/5 + i];\
		m01234[i] = m0123 + matr[8*X/5 + i];\
		m50123[i] = matr[3*X/5 + i] + m0123;\
		m65 = matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = matr[1*X/5 + i] + m65;\
		m65012[i] = m65 + m012;\
		m76501[i] = m765 + m01;\
		m87650[i] = matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	F (t5 , v0m4, matr + 4*X/5);\
	F (t6 , v1m4, matr + 3*X/5);\
	F (t7 , v2m4, matr + 2*X/5);\
	F (t8 , v3m4, matr + 1*X/5);\
	F (t9 , v0m3, matr + 5*X/5);\
	F (t10, v1m3, matr + 4*X/5);\
	F (t11, v2m3, matr + 3*X/5);\
	F (t12, v0m2, matr + 6*X/5);\
	F (t13, v1m2, matr + 5*X/5);\
	F (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#define bigmatTOEP55TOP(X, F) {\
	__int128 t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	int64_t v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	__int128 m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	__int128 m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = vect[i] - vect[4*X/5 + i];\
		v1m4[i] = vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = vect[i] - vect[3*X/5 + i];\
		v1m3[i] = vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = vect[i] - vect[2*X/5 + i];\
		v1m2[i] = vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = (__int128)matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = (__int128)m01 + matr[6*X/5 + i];\
		m0123 = (__int128)m012 + matr[7*X/5 + i];\
		m01234[i] = (__int128)m0123 + matr[8*X/5 + i];\
		m50123[i] = (__int128)matr[3*X/5 + i] + m0123;\
		m65 = (__int128)matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = (__int128)matr[1*X/5 + i] + m65;\
		m65012[i] = (__int128)m65 + m012;\
		m76501[i] = (__int128)m765 + m01;\
		m87650[i] = (__int128)matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	F (t5 , v0m4, matr + 4*X/5);\
	F (t6 , v1m4, matr + 3*X/5);\
	F (t7 , v2m4, matr + 2*X/5);\
	F (t8 , v3m4, matr + 1*X/5);\
	F (t9 , v0m3, matr + 5*X/5);\
	F (t10, v1m3, matr + 4*X/5);\
	F (t11, v2m3, matr + 3*X/5);\
	F (t12, v0m2, matr + 6*X/5);\
	F (t13, v1m2, matr + 5*X/5);\
	F (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#define M1SCHOOLBOOK(X) {\
	for(int i = 0; i < X; i++)\
	{\
		rop[i] = 0;\
		for(int j = 0; j < X; j++)\
			rop[i] += vect[j] * matr[X - 1 - j + i];\
	}\
}

#define TOEP55TOP128(X, F) {\
	__int128 t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	__int128 v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	__int128 m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	__int128 m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = (__int128)vect[i] - vect[4*X/5 + i];\
		v1m4[i] = (__int128)vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = (__int128)vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = (__int128)vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = (__int128)vect[i] - vect[3*X/5 + i];\
		v1m3[i] = (__int128)vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = (__int128)vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = (__int128)vect[i] - vect[2*X/5 + i];\
		v1m2[i] = (__int128)vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = (__int128)vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = (__int128)matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = (__int128)m01 + matr[6*X/5 + i];\
		m0123 = (__int128)m012 + matr[7*X/5 + i];\
		m01234[i] = (__int128)m0123 + matr[8*X/5 + i];\
		m50123[i] = (__int128)matr[3*X/5 + i] + m0123;\
		m65 = (__int128)matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = (__int128)matr[1*X/5 + i] + m65;\
		m65012[i] = (__int128)m65 + m012;\
		m76501[i] = (__int128)m765 + m01;\
		m87650[i] = (__int128)matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	F (t5 , v0m4, matr + 4*X/5);\
	F (t6 , v1m4, matr + 3*X/5);\
	F (t7 , v2m4, matr + 2*X/5);\
	F (t8 , v3m4, matr + 1*X/5);\
	F (t9 , v0m3, matr + 5*X/5);\
	F (t10, v1m3, matr + 4*X/5);\
	F (t11, v2m3, matr + 3*X/5);\
	F (t12, v0m2, matr + 6*X/5);\
	F (t13, v1m2, matr + 5*X/5);\
	F (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#define M1SCHOOLBOOK(X) {\
	for(int i = 0; i < X; i++)\
	{\
		rop[i] = 0;\
		for(int j = 0; j < X; j++)\
			rop[i] += vect[j] * matr[X - 1 - j + i];\
	}\
}

#define M1TOEP22TOP(X, F) {\
	int64_t t0[X/2], t1[X/2], t2[X/2];\
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

#define M1TOEP33TOP(X, F) {\
	int64_t t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3];\
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

#define M1TOEP55TOP(X, F) {\
	int64_t t0[X/5], t1[X/5], t2[X/5], t3[X/5], t4[X/5], t5[X/5], t6[X/5], t7[X/5], t8[X/5], t9[X/5], t10[X/5], t11[X/5], t12[X/5], t13[X/5], t14[X/5];\
	int64_t v0m4[X/5], v1m4[X/5], v2m4[X/5], v3m4[X/5], v0m3[X/5], v1m3[X/5], v2m3[X/5], v0m2[X/5], v1m2[X/5], v0m1[X/5];\
	int64_t m01234[(2*X/5) - 1], m50123[(2*X/5) - 1], m65012[(2*X/5) - 1], m76501[(2*X/5) - 1], m87650[(2*X/5) - 1];\
	int64_t m01, m012, m0123, m65, m765;\
	for(int i = 0; i < X/5; i++)\
	{\
		v0m4[i] = vect[i] - vect[4*X/5 + i];\
		v1m4[i] = vect[X/5 + i] - vect[4*X/5 + i];\
		v2m4[i] = vect[2*X/5 + i] - vect[4*X/5 + i];\
		v3m4[i] = vect[3*X/5 + i] - vect[4*X/5 + i];\
		v0m3[i] = vect[i] - vect[3*X/5 + i];\
		v1m3[i] = vect[X/5 + i] - vect[3*X/5 + i];\
		v2m3[i] = vect[2*X/5 + i] - vect[3*X/5 + i];\
		v0m2[i] = vect[i] - vect[2*X/5 + i];\
		v1m2[i] = vect[X/5 + i] - vect[2*X/5 + i];\
		v0m1[i] = vect[i] - vect[X/5 + i];\
	}\
	for(int i = 0; i < (2*X/5) - 1; i++)\
	{\
		m01 = matr[4*X/5 + i] + matr[5*X/5 + i];\
		m012 = m01 + matr[6*X/5 + i];\
		m0123 = m012 + matr[7*X/5 + i];\
		m01234[i] = m0123 + matr[8*X/5 + i];\
		m50123[i] = matr[3*X/5 + i] + m0123;\
		m65 = matr[2*X/5 + i] + matr[3*X/5 + i];\
		m765 = matr[1*X/5 + i] + m65;\
		m65012[i] = m65 + m012;\
		m76501[i] = m765 + m01;\
		m87650[i] = matr[i] + m765 + matr[4*X/5 + i];\
	}\
	F (t0, vect + 4*X/5, m87650);\
	F (t1, vect + 3*X/5, m76501);\
	F (t2, vect + 2*X/5, m65012);\
	F (t3, vect + 1*X/5, m50123);\
	F (t4, vect        , m01234);\
	F (t5 , v0m4, matr + 4*X/5);\
	F (t6 , v1m4, matr + 3*X/5);\
	F (t7 , v2m4, matr + 2*X/5);\
	F (t8 , v3m4, matr + 1*X/5);\
	F (t9 , v0m3, matr + 5*X/5);\
	F (t10, v1m3, matr + 4*X/5);\
	F (t11, v2m3, matr + 3*X/5);\
	F (t12, v0m2, matr + 6*X/5);\
	F (t13, v1m2, matr + 5*X/5);\
	F (t14, v0m1, matr + 7*X/5);\
	for(int i = 0; i < X/5; i++)\
	{\
		rop[i]         = t0[i] + t5[i] + t6[i]  + t7[i]  + t8[i];\
		rop[i +   X/5] = t1[i] - t8[i] + t9[i]  + t10[i] + t11[i];\
		rop[i + 2*X/5] = t2[i] - t7[i] - t11[i] + t12[i] + t13[i];\
		rop[i + 3*X/5] = t3[i] - t6[i] - t10[i] - t13[i] + t14[i];\
		rop[i + 4*X/5] = t4[i] - t5[i] - t9[i]  - t12[i] - t14[i];\
	}\
}

#endif

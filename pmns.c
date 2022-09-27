#include <stdio.h>
#include <time.h>

#include "pmns.h"
#include "utilitymp.h"

#define AMNS_MONTG_MULT amns_montg_mult_pre

__inline uint64_t mulx64(uint64_t x, uint64_t y, uint64_t* hi)
{
    __asm__(
        "mulx %3, %0, %1    \n\t"
        : "=&d"(x), "=&a"(y)
        : "0"(x), "1"(y)
    );

    *hi = y;
    return x;
}

inline void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	
/*	uint64_t *Rhi, *Rlo;*/
/*	__int128 Rfull;*/
/*	*/
/*	Rlo = (uint64_t*) &Rfull;*/
/*	Rhi = Rlo + 1;*/
/*	uint64_t hi, lo;*/
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
/*		{*/
/*			*Rlo = mulx64(A->t[i + j], B->t[N - j], Rhi);*/
/*			R[i] += Rfull;*/
/*			lo = mulx64(A->t[i + j], B->t[N - j], &hi);*/
/*			R[i] += (__int128) lo | ((__int128) hi << 64);*/
/*		}*/
		
		R[i] = R[i] * LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
/*		{*/
/*			*Rlo = mulx64(A->t[j], B->t[i - j], Rhi);*/
/*			R[i] += Rfull;*/
/*			lo = mulx64(A->t[j], B->t[i - j], &hi);*/
/*			R[i] += (__int128) lo | ((__int128) hi << 64);*/
/*		}*/
	}
}

static inline void m_mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) (A->t[i + j]) * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) (A->t[j]) * M[i - j];
	}
}

inline void m1_mns_mod_mult_ext_red(int64_t* restrict R,
	const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += ((uint64_t)A->t[i + j]) * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += ((uint64_t)A->t[j]) * M1[i - j];
	}
}

inline void mns_montg_int_red(restrict poly res, __int128* restrict R)
{
	uint64_t V[N], V2[N];
	int64_t T[N] = {0};
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	m1_mns_mod_mult_ext_red(T, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = T[i];
	
	m_mns_mod_mult_ext_red(R, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = V2[i] + (R[i] >> 64) + (V[i] != 0);
}

inline void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 R[N] = {0};
	
	mns_mod_mult_ext_red(R, A, B);

	mns_montg_int_red(res, R);
}

void amns_rtl_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		for(j = 0; j < 64; j++)
			{
				if(aux & (1ULL << j))
					AMNS_MONTG_MULT(tmp, res, base);
				else
					AMNS_MONTG_MULT(tmp, res, &__theta__);
				AMNS_MONTG_MULT(res, tmp, tmp);
			}
	}
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			AMNS_MONTG_MULT(tmp, res, base);
		else
			AMNS_MONTG_MULT(tmp, res, &__theta__);
		AMNS_MONTG_MULT(res, tmp, tmp);
		aux >>= 1;
	}
	
	free_poly(tmp);
}

void amns_ltr_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux) + 1; j < 64; j++)
	{
		AMNS_MONTG_MULT(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			AMNS_MONTG_MULT(res, tmp, base);
		else
			AMNS_MONTG_MULT(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			AMNS_MONTG_MULT(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				AMNS_MONTG_MULT(res, tmp, base);
			else
				AMNS_MONTG_MULT(res, tmp, &__theta__);
		}
	}
	
	free_poly(tmp);
}

void amns_sqandmult(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	amns_rtl_sqandmult(res, base, exponent);
}

void amns_montg_ladder(restrict poly res, const restrict poly base,
	const restrict poly exponent)
{
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux, b;
	
	poly tmp, R[2];
	init_poly(N, &tmp);
	R[0] = res;
	R[1] = tmp;
	
	poly_copy(res, &__theta__);
	poly_copy(tmp, base);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux) + 1; j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		AMNS_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
		AMNS_MONTG_MULT(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			AMNS_MONTG_MULT(R[1 - b], R[1 - b], R[b]);
			AMNS_MONTG_MULT(R[b], R[b], R[b]);
		}
	}
	
	free_poly(tmp);
}

static inline void mns_montg_int_red_pre(restrict poly res, __int128* restrict R)
{
	uint64_t V[N], V2[N];
	int64_t T[N] = {0};
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	m1_mns_mod_mult_ext_red_pre(T, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = T[i];
	
	m_mns_mod_mult_ext_red_pre(R, res);
	
	for(i = 0; i < N; i++)
		res->t[i] = V2[i] + (R[i] >> 64) + (V[i] != 0);
}

inline void amns_montg_mult_pre(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	
	__int128 R[N] = {0};
	
	mns_mod_mult_ext_red_pre(R, A, B);

	mns_montg_int_red_pre(res, R);
}

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

static inline int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

void randpoly(poly P)
{
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64()) * (1 + (rand() & 1) * -2);
}

void convert_string_to_amns(restrict poly res, const char* string)
{
	uint8_t counter;
	register uint16_t i, j;
	const uint64_t rho = (1ULL<<RHO);
	__int128 R[N] = {0};
	poly stok;
	init_poly(N, &stok);
	
	convert_string_to_poly(&stok, string);
	
	if(stok->deg > N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	for(i = 1; i < N; i++)
		res->t[i] = 0;
	
	res->t[0] = ((uint64_t) stok->t[0]) & (rho - 1);
	counter = 0;
	for(i = 1; i < N; i++)
	{
		counter = (counter + 64 - RHO);
		res->t[i] += ((((uint64_t) stok->t[i]) << counter) | (((uint64_t) stok->t[i - 1]) >> (64 - counter))) & (rho - 1);
		if(counter > RHO && i < N - 1)
		{
			counter = counter - RHO;
			res->t[i + 1] = ((((uint64_t) stok->t[i]) << counter) | (((uint64_t) stok->t[i - 1]) >> (64 - counter))) & (rho - 1);
		}
	}
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);

end:
	free_poly(stok);
}

void convert_amns_to_poly(restrict poly* res, const restrict poly P)
{
	register uint16_t i;
	poly a, aux, ag, tmp;
	__int128 Quite[N];
	
	init_polys(N, &a, &ag, &tmp, NULL);
	init_poly(1, &aux);
	if((*res)->deg < N)
	{
		free_poly(*res);
		init_poly(N, res);
	}
	(*res)->deg = N;
	for(i = 1; i < N; i++)
	{
		(*res)->t[i] = 0;
		Quite[i] = (__int128) P->t[i];
	}
	Quite[0] = (__int128) P->t[0];
	
	mns_montg_int_red(a, Quite);
	
	(*res)->t[0] = a->t[0];
	for(i = 1; i < N; i++)
	{
		aux->t[0] = a->t[i];
		mp_mult(&ag, aux, &Gi[i - 1]);
		
		mp_copy(&tmp, *res);
		mp_add(res, tmp, ag);
	}
	
	mp_copy(&tmp, *res);
	mp_mod(res, tmp, &__P__);
	
	free_polys(a, aux, ag, tmp, NULL);
}

void __sub_tests(void)
{
	poly A, B, AmB, BmA, AmBmBmA, BmAmAmB, ApB, AmBpBmA, BmApAmB;
	init_polys(5, &A, &B, &AmB, &BmA, &AmBmBmA, &BmAmAmB, NULL);
	init_polys(6, &ApB, &AmBpBmA, &BmApAmB, NULL);

	const char a[] = "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
		b[] = "c64be1b07fc889170737a3bbd0501940eb5cdffbb228d09f6c4527812aa64b8cde15";
	convert_string_to_poly(&A, a);
	mp_print(A);
	convert_string_to_poly(&B, b);
	mp_print(B);
	mp_sub(&AmB, A, B);
	mp_print(AmB);
	mp_sub(&BmA, B, A);
	mp_print(BmA);
	mp_sub(&AmBmBmA, AmB, BmA);
	mp_print(AmBmBmA);
	mp_sub(&BmAmAmB, BmA, AmB);
	mp_print(BmAmAmB);
	mp_add(&ApB, A, B);
	mp_print(ApB);
	mp_add(&BmApAmB, BmA, AmB);
	mp_print(BmApAmB);
	mp_add(&AmBpBmA, AmB, BmA);
	mp_print(AmBpBmA);
	printf("%d %d %d %d\n", mp_comp(A, B), mp_comp(A, AmB), mp_comp(ApB, B),
		mp_comp(AmBpBmA, BmApAmB));
	printf("-1 1 1 0\n");

	free_polys(A, B, AmB, BmA, AmBmBmA, BmAmAmB, ApB, AmBpBmA, BmApAmB, NULL);
}

void __mult_tests(void)
{
	poly A, B, C, C_test;
	init_polys(1, &A, &B, NULL);
	init_polys(1, &C, &C_test, NULL);
	
	const char a[] =  "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
		b[] =  "c64be1b07fc889170737a3bbd0501940eb5cdffbb228d09f6c4527812aa64b8cde15",
		c[] =  "113a594a40463202213fb6a4f7f7c643b20e4efcf7206e5154247afe315d4d6c9a5f8dcd34573cc1b63874bccb81c44e190796da9afb9b7f144f1fe35456d916cebebdf9";
	
	convert_string_to_poly(&A, a);
	printf("%s\n", a);
	mp_print(A);
	convert_string_to_poly(&B, b);
	printf("%s\n", b);
	mp_print(B);
	convert_string_to_poly(&C_test, c);
	printf("%s\n", c);
	mp_print(C_test);
	mp_mult(&C, A, B);
	mp_print(C);
		
	free_polys(A, B, C, C_test, NULL);
}

void __mod_tests(void)
{
	poly A, B, C, ST;
	init_polys(0, &A, &B, &C, &ST, NULL);
	
	const char a[] = "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
	b[] = "a84c8a29612d9545394797803178257e7e72af34d038e704ca1d0e2000ad0f01e8e03f441888cdbcf393542b3f896520e3f53e37e98e4e7dfe83645e74d239315d054901",
	c[] = "be42745e4e6c627ae20e738fdf363efa5ad0c8cc844644528d71ac55ab9d69baa36";
	convert_string_to_poly(&A, a);
	mp_print(A);
	convert_string_to_poly(&B, b);
	mp_print(B);
	convert_string_to_poly(&ST, "1");
	for(int i = 0; i < 75; i++)
	{
		mp_leftshift(&ST);
		ST->t[0] |= 1;
	}
	for(int i = 0; i < 75; i++)
		mp_rightshift(ST);
	
	mp_mod(&C, B, A);
	printf("%s\n", c);
	mp_print(C);
	
	free_polys(A, B, C, ST, NULL);
}

void __init_tests__(void)
{
	poly a, b, c, c_check, Phisquared;
	init_polys(N, &a, &b, &c, &c_check, &Phisquared, NULL);
	
	const char A[] = "9b6afe91a6e17ff3e5b7331bc220a825e6bbe48687ca568a0873800b48471d633375";
	poly converted;
	init_poly(N, &converted);
	convert_string_to_amns(converted, A);
	printf("%s\n", A);
	print(converted);
	free_poly(converted);	
	
	set_val(a, 3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128);
	set_val(b, 1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068);
	//set_val(b, 1, 0, 0, 0, 0);
	//set_val(c_check, 2302327877203981, 25683149970777821, -1798382075251775, 52479742770215631, 21994577573493812);
	set_val(Phisquared, 0, 0, 0, 512, 0);
	
	//amns_montg_mult(c, a, b);
	amns_montg_mult(c_check, a, Phisquared);
	amns_montg_mult(c, c_check, b);
	
	//print(Mlambda);
	//print(M1lambda);
	print(a);
	print(b);
	//print(c_check);
	print(c);
	
	free_polys(a, b, c, c_check, Phisquared, NULL);
}

void __proof_of_accuracy(void)
{
	time_t seed;
	poly a, aphi, b, c, Phisquared;
	init_polys(N, &a, &aphi, &b, &c, &Phisquared, NULL);
	set_val(Phisquared, 0, 0, 0, 512, 0);
	srand((unsigned)time(&seed));
	
		for(register int64_t i = 0; i < 100; i++)
	{
		randpoly(a);
		amns_montg_mult(aphi, a, Phisquared);
		randpoly(b);
		amns_montg_mult(c, aphi, b);
		
		print(a);
		print(b);
		print(c);
	}
	
	free_polys(a, aphi, b, c, Phisquared, NULL);
}

void __full_mult_demo(void)
{
	const char a[] = "77f882926258fb5a293015e16fc961598939f9f328d4e316d02519d3f8d88412d787",
		b[] = "b4399ccbab87f4f053d75a9dcc1c1fa8d2f4edd7bdf5eebc78fb4ea16a6fb02eb96d",
		c[] = "5475bb9ee1c1ea9bb15ef6fc9118d391a012ab14a4ae063bc8832426a1e5c0f555a8fb1e335846abcf0e969b6c942fe0bf34126cfb7f1477a28dc2147f1a74ba6408537b";
	
	poly A, B, C, aux;
	
	init_polys(N, &A, &B, &C, &aux, NULL);
	
	convert_string_to_amns(A, a);
	convert_string_to_amns(B, b);
	amns_montg_mult(C, A, B);
	
	convert_amns_to_poly(&aux, C);
	
	poly tmp;
	init_poly(1, &tmp);
	
	_poly tmp1 = { .deg = 1, .t = (int64_t[]) { -10 } },
		tmp2 = { .deg = 1, .t = (int64_t[]) { 7 } };
	
	mp_mod(&tmp, &tmp1, &tmp2);
	
	mp_print(tmp);
	free_poly(tmp);
	
	printf("\n\n");
	
	mp_print(aux);
	printf("%s\n", c);
	
	free_polys(A, B, C, aux, NULL);
}

void __multbench__(void)
{
	uint64_t c1 = 0, sum;
	poly a, b, c, soak1, soak2;
	init_polys(N, &a, &b, &c, &soak1, &soak2, NULL);
	
	srand((unsigned) (time(((int64_t*)(&c1)))));
	
	c1 = clock();
	
	randpoly(soak2);
	soak2->t[0] += Gi[0].t[0];
	soak2->t[0] += __P__.t[0];
	
	for(int i = 0; i < 1000; i++)
	{
		randpoly(a);
		randpoly(b);
		amns_montg_mult(c, a, b);
		amns_montg_mult(soak1, c, soak2);
		amns_montg_mult(soak2, c, soak1);
	}
	
	sum = 0;
	for(int i = 0; i < 10000; i++)
	{
		randpoly(a);
		randpoly(b);
		c1 = clock();
		amns_montg_mult(c, a, b);
		sum += clock() - c1;
		amns_montg_mult(soak1, c, soak2);
		amns_montg_mult(soak2, c, soak1);
	}
	
	printf("%ld\n", sum);
	
	free_polys(a, b, c, soak1, soak2, NULL);
}

void __multchecks__(void)
{
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	int64_t seed;
	
	srand((unsigned) (time(&seed)));
	
	for(int i = 0; i < 100; i++)
	{
		randpoly(a);
		randpoly(b);
		print(a);
		print(b);
		amns_montg_mult(c, a, b);
		print(c);
	}
	
	free_polys(a, b, c, NULL);
}

void __sqandmultdemo(void)
{
	const char a[] = "2", // "7609d69beaadb6de37a7b36cd193b33b120489bd4298534e830eeeaf9a65b15c12268aea1447f610377ea045afc463fb193a531e46cf70052ee6143d782b27aee363d426ad73085f7c24376b676070214cbf1b69f93fd5fdd70b8c77dd2268cbf3f210366b932c7351d9332608fb294ebb44bc7b17bfa3e115dd06c642670d67",
		b[] = "1", // "78a5cc1a942f42e81aa0dc980d3b6a7f987bf9847fb30f8d17c327283745d1365e6bd495a6b4bc2f6a16dca99668ee8591b5b04d12e15d5c4f5f22aafca94fcc3973e7e7714c84b0fb514f862d3444fd36f82f54ccb6c2f2ac3510d3aaea94953533e511076ba29103afb13ba6387001f1055b400b5abfc5b7cd0cdb4b2ebf2a",
		c[] = "606ddc8f63ff63abfacb0ced7fcef39d42f0831e9e84459bbffbc1a04b866aaf572571e876a087dc633ab189b40eb861be2f7e194ac5f24ad7886eeb070028bc91a3970ea828cdd5da8eea173d0a38da1bc072837e5835ff2bd266ebcc4788870c7dca82e5b87cd1a844a5120339d2ef1620f1484888538824c5474fe77014e6";
	
	poly A, B, C, aux;
	init_polys(N, &A, &B, &C, &aux, NULL);
	
	convert_string_to_poly(&B, b);
	convert_string_to_amns(A, a);
	
	//amns_sqandmult(C, A, B);
	amns_montg_ladder(C, A, B);
	
	//printf("C = ");
	//print(C);
	
	convert_amns_to_poly(&aux, C);
	
	printf("\n%s\n\n", c);
	mp_print(aux);
	
	free_polys(A, B, C, aux, NULL);
}

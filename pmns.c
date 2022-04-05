#include <stdio.h>
#include <time.h>

#include "pmns.h"
#include "utilitymp.h"



static inline void mns_mod_mult_ext_red(__int128* R, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * B->t[N - j];
		
		R[i] = R[i] * LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

static inline void m_mns_mod_mult_ext_red(__int128* R, const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) ((uint64_t)A->t[i + j]) * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) ((uint64_t)A->t[j]) * M[i - j];
	}
}

static inline void m1_mns_mod_mult_ext_red(__int128* R, const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M1 being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) ((uint64_t)A->t[i + j]) * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) ((uint64_t)A->t[j]) * M1[i - j];
	}
}

inline void mns_montg_int_red(restrict poly res, __int128* R)
{
	uint64_t V[N], V2[N], T[N], T2[N];
	register uint16_t i;
	
	for(i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	
	m1_mns_mod_mult_ext_red(R, res);
	
	for(i = 0; i < N; i++)
	{
		res->t[i] = R[i];
		R[i] = 0;
	}
	
	m_mns_mod_mult_ext_red(R, res);
	
	for(i = 0; i < N; i++)
	{
		T[i] = R[i];
		T2[i] = (R[i] >> 64);
		
		T[i] += V[i];
		res->t[i] = V2[i] + T2[i] + (T[i] < V[i]);
	}
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

static inline int64_t randomint64(void)
{
	return (int64_t)(((int64_t)(rand() + rand()) << 32) ^ ((int64_t)(rand() + rand())));
}

static inline int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

static inline void randpoly(poly P)
{
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64());
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
	init_polys(5, &A, &B, NULL);
	init_polys(10, &C, &C_test, NULL);
	
	const char a[] = "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
		b[] = "c64be1b07fc889170737a3bbd0501940eb5cdffbb228d09f6c4527812aa64b8cde15",
		c[] = "113a594a40463202213fb6a4f7f7c643b20e4efcf7206e5154247afe315d4d6c9a5f8dcd34573cc1b63874bccb81c44e190796da9afb9b7f144f1fe35456d916cebebdf9";
	
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
	
		for(register int64_t i = 0; i < 100000; i++)
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
/*		c[] = "199da60a1e7f5c5523c51230e885506efe15b08b191b981590b47c0567b3bb516da3";*/
		c[] = "129cc51d3f86bb90b6cdf86018bf9bbb6353602a99c70e7e90305551f8e9d88c288277746a878ec690bad018212c4e7eaa2577215916addc79b7842c303898d5cf6";
	
	poly A, B, C, aux;
	
	init_polys(N, &A, &B, &C, &aux, NULL);
	
	convert_string_to_amns(A, a);
	convert_string_to_amns(B, b);
	amns_montg_mult(C, A, B);
	
	convert_amns_to_poly(&aux, C);
	
	mp_print(aux);
	printf("%s\n", c);
	
	free_polys(A, B, C, aux, NULL);
}


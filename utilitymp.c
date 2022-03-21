#include <stdio.h>
#include "montgom.h"
#include "utilitymp_core.h"


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


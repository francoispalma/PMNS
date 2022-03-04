#include <stdio.h>
#include <stdint.h>

#include "montgom.h"

static inline void printasonenumber(poly P)
{
	printf("%lx", P->t[P->deg - 1]);
	for(register uint16_t i = 1; i < P->deg; i++)
		printf("%016lx", P->t[P->deg - 1 - i]);
	printf("\n"); 
}

static inline void multmp(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	register int16_t i, j, k;
	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < op1->deg + op2->deg)
	{
		free_poly(*res);
		init_poly(op1->deg + op2->deg, res);
	}
	
	__int128 R[op1->deg + op2->deg], stok;
	for(i = 0; i < op1->deg + op2->deg; i++)
		R[i] = 0;
	
	for(i = 0; i < op1->deg; i++)
	{
		for(j = 0; j < op2->deg; j++)
		{
			k = i + j;
			stok = R[k];
			R[k] += (__int128) op1->t[i] * op2->t[j];
			while(stok > R[k])
			{
				++k;
				stok = R[k];
				R[k] += 1;
			}
		}
	}
	
	res->t[0] = R[0];
	for(i = 1; i < op1->deg + op2->deg; i++)
	{
		k = i;
		stok = R[i];
		R[i] += R[i - 1] >> 64;
		while(stok > R[k])
		{
			++k;
			stok = R[k];
			R[k] += 1;
		}
		res->t[i] = R[i];
	}
}

int main(void)
{
	poly A, B, C, C_test;
	init_polys(5, &A, &B, NULL);
	init_polys(10, &C, &C_test, NULL);
	
	const char a[] = "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
		b[] = "c64be1b07fc889170737a3bbd0501940eb5cdffbb228d09f6c4527812aa64b8cde15",
		c[] = "113a594a40463202213fb6a4f7f7c643b20e4efcf7206e5154247afe315d4d6c9a5f8dcd34573cc1b63874bccb81c44e190796da9afb9b7f144f1fe35456d916cebebdf9";
	
	convert_string_to_poly(&A, a);
	printf("%s\n", a);
	printasonenumber(A);
	convert_string_to_poly(&B, b);
	printf("%s\n", b);
	printasonenumber(B);
	convert_string_to_poly(&C_test, c);
	printf("%s\n", c);
	printasonenumber(C_test);
	multmp(C, A, B);
	printasonenumber(C);
		
	free_polys(A, B, C, C_test, NULL);
	return 0;
}

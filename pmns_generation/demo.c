#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#include "structs.h"
#include "params.h"
#include "pmns.h"
#include "utilitymp.h"

#define UNUSED(X) ((void)X)

static inline int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

int main(void)
{
	int64_t seed;
	
	srand((unsigned) (time(&seed)));
	
	mpnum A; poly AA; mpnum AAA; mpnum tmp;
	init_mpnums(N, &A, &AAA, &tmp, 0);
	init_poly(N, &AA);
	for(int i = 0; i < N; i++)
		tmp->t[i] = randomint64();
	mp_mod(&A, tmp, (const mpnum) &__P__);
	printf("Random element from Z/pZ:\n");
	mp_print(A);
	printf("\nPolynomial conversion to our PMNS:\n");
	convert_binary_to_pmns(AA, A);
	poly_print(AA);
	printf("\nConverting back to binary:\n");
	convert_pmns_to_binary(&AAA, AA);
	mp_print(AAA);
	
	printf("\nEquality between the two:\n");
	if(mp_comp(A, AAA) == 0)
		printf("True\n");
	else
		printf("False\n");
	
	free_mpnums(A, AAA, tmp, 0);
	free_poly(AA);
	
	UNUSED(Gi);
	UNUSED(__theta__);
	
	return 0;
}

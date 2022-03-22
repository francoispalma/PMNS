#include <stdio.h>

#include "mppmns.h"
#include "utilitymp.h"

#define LOW(X) ((uint64_t)X)
#define HIGH(X) ((uint64_t)(X>>64))

static inline void mns128_montg_int_red(poly128 res, const __int128* R)
{
	// pass
}

void convert_string_to_amns128(restrict poly128 res, const char* string)
{
	uint8_t counter;
	register uint16_t i, j;
	const unsigned __int128 rho = ((__int128)1) << RHO;
	unsigned __int128 limb;
	__int128 R[N] = {0}, tmp[N] = {0};
	poly stok;
	init_poly(2 * N, &stok);
	
	convert_string_to_poly(&stok, string);
	printf("%s\n", string);
	mp_print(stok);
	
	if(stok->deg > 2 * N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	counter = 0;
	for(i = 0; i < N - 1; i++)
	{
		limb = (((unsigned __int128) stok->t[2 * i + 1]) << 64)
		 | ((uint64_t) stok->t[2 * i]);
		tmp[i] |= (limb << counter) & (rho - 1);
		__print128(tmp[i]);
		tmp[i + 1] |= (limb >> (RHO - counter));
		counter = (counter + 128 - RHO) % RHO;
	}
	limb = (unsigned __int128) ((unsigned __int128) stok->t[2 * i]) |
			(((unsigned __int128) stok->t[2 * i + 1]) << 64);
	tmp[i] |= (limb << counter) & (rho - 1);
	__print128(tmp[i]);
	
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) tmp[i] * __Pilo__[i][j];
	
	mns128_montg_int_red(res, R);
	
end:
	free_poly(stok);
}


int main(void)
{
	const char a[] = "74ff560400d0105e6381e4f7cf22ba4a3d949bbe3b03e7ec1c8aebfb02a4dedf230eef099cd1ae78adf8f142cd70ed93122a5c48c5edcba658615fa2316994dce0c84e9e54c5ae9482acdc0ed6fae84eb7e83d94016d12452ad41369e33a53a676d539439488bdc8b3462c5579a432e8b579e8af9d5b2b0b8f37856fe2de7f30";
	
	convert_string_to_amns128(NULL, a);
	
	return 0;
}


/*
['0x39e8af9d5b2b0b8f37856fe2de7f30', '0x250e5222f722cd18b155e690cba2d5', '0x16d12452ad41369e33a53a676d53', '0x316ba520ab3703b5beba13adfa0f65', '0xba658615fa2316994dce0c84e9e54', '0x22b7e3c50b35c3b64c48a9712317b7', '0x8aebfb02a4dedf230eef099cd1ae7', '0x393df3c8ae928f6526ef8ec0f9fb07', '0x74ff560400d0105e6381']
*/

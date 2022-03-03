#include <stdio.h>
#include <time.h>

#include "montgom.h"

#define RHO 61
#define N 5
#define LAMBDA 2

// Various values you can get from precalc
const int64_t M[] = {9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940}, M1[] = {7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353}, MLambda[] = {18892189295192050, -178689719550530, 732756003058628, -9116351660287002, 38462084564469880}, M1Lambda[] = {-4355481239625866354, 6277017080637277974, -2433945398614454754, 3202559735019019384, 8710962479251732706}, __Pi__[N][N] = {
	{13165649034348498, 26180244800938503, 13254993894123763, 25813866799409701, 17813169724267264},
	{27785276313683301, 7751718275075962, 5780308760799802, 15841760209993152, 16644275966365551},
	{48640935215124533, 240398127541950, 38739117706223929, 5254296818981838, 14583521592316174}, 
	{31683520419953536, 33288551932731102, 27785276313683301, 7751718275075962, 5780308760801850}, 
	{13165649034348498, 26180244800938503, 13254993894123763, 25813866799409189, 17813169724332800}
};



inline void init_poly(const uint16_t deg, restrict poly* P)
{
	*P = malloc(sizeof(_poly));
	(*P)->deg = deg;
	(*P)->t = calloc(deg, sizeof(int64_t));
}

void init_polys(const uint16_t deg, restrict poly* P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		init_poly(deg, P);
		P = va_arg(args, poly*);
	} while(P != NULL);
	va_end(args);
}

inline void free_poly(restrict poly P)
{
	free(P->t);
	free(P);
}

void free_polys(restrict poly P, ...)
{
	va_list args;
	
	va_start(args, P);
	do
	{
		free_poly(P);
		P = va_arg(args, poly);
	} while(P != NULL);
	va_end(args);
}

void set_val(restrict poly P, int64_t val, ...)
{
	va_list args;
	
	va_start(args, val);
	for(int16_t i = 0; i < P->deg; i++)
	{
		P->t[i] = val;
		val = va_arg(args, int64_t);
	}
	va_end(args);
}

inline void print(const restrict poly P)
{
	printf("[");
	for(int16_t i = 0; i < P->deg - 1; i++)
		printf("%ld, ", P->t[i]);
	printf("%ld]\n", P->t[P->deg - 1]);
}

static inline void __print128(register const __int128 Val)
{
	int64_t hi = Val >> 64;
	uint64_t lo = Val;
	if (hi & 0x8000000000000000) printf("-");
	printf("0x%lx%016lx\n", hi & 0x8000000000000000 ? -hi : hi, lo);
}

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
			R[i] += (__int128) A->t[i + j] * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * M[i - j];
	}
}

static inline void m1_mns_mod_mult_ext_red(__int128* R, const restrict poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A->t[i + j] * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * M1[i - j];
	}
}

static inline void mns_montg_int_red(restrict poly res, __int128* R)
{
	uint64_t V[N], V2[N], T[N], T2[N];

	for(int i = 0; i < N; i++)
	{
		V[i] = R[i];
		res->t[i] = R[i];
		V2[i] = (R[i] >> 64);
		R[i] = 0;
	}
	m1_mns_mod_mult_ext_red(R, res);
	
	for(int i = 0; i < N; i++)
	{
		res->t[i] = R[i];
		R[i] = 0;
	}
	
	m_mns_mod_mult_ext_red(R, res);
	
	for(int i = 0; i < N; i++)
	{
		T[i] = R[i];
		T2[i] = (R[i] >> 64);
		
		T[i] = V[i] + T[i];
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

inline void convert_string_to_poly(restrict poly* res, const char* string)
{
	uint16_t tabsize = 0, j;
	int16_t i = 0;
	char store[17];
	
	store[16] = '\0';
	
	while(string[i] != '\0')
	{
		if(i % 16 == 0) ++tabsize;
		++i;
	}
	if (tabsize > (*res)->deg)
	{
		free_poly(*res);
		init_poly(tabsize, res);
	}

	// we keep i from the end of the last loop on purpose.
	for(i=i; i > -1; i -= 16)
	{
		for(j = 0; j < 16; j++)
		{
			if(i - 16 + j >= 0)
				store[j] = string[i - 16 + j];
			else
				store[j] = '0';
		}
		(*res)->t[tabsize - 1 - (i / 16)] = strtoul(store, NULL, 16);
	}
}

void convert_string_to_amns(restrict poly res, const char* string)
{
/*	uint16_t tabsize = 0, maxi, j;*/
/*	int16_t i = 0;*/
/*	uint8_t counter;*/
/*	const uint64_t rho = (1ULL<<RHO);*/
/*	uint64_t tab[N];*/
/*	__int128 R[N] = {0};*/
/*	char store[17];*/
/*	*/
/*	store[16] = '\0';*/
/*	*/
/*	while(string[i] != '\0')*/
/*	{*/
/*		if(i % 16 == 0) ++tabsize;*/
/*		++i;*/
/*	}*/
/*	maxi = i;*/
/*	for(i = maxi; i > -1; i -= 16)*/
/*	{*/
/*		for(j = 0; j < 16; j++)*/
/*		{*/
/*			if(i - 16 + j >= 0)*/
/*				store[j] = string[i - 16 + j];*/
/*			else*/
/*				store[j] = '0';*/
/*		}*/
/*		tab[tabsize - 1 - (i / 16)] = strtoul(store, NULL, 16);*/
/*	}*/
	


	uint8_t counter;
	uint16_t i, j;
	const uint64_t rho = (1ULL<<RHO);
	__int128 R[N];
	poly stok;
	init_poly(N, &stok);
	
	convert_string_to_poly(&stok, string);
	printf("%lx", stok->t[N - 1]);
	for(i = 1; i < N; i++)
		printf("%016lx", stok->t[N - 1 - i]);
	printf("\n");
	//TODO check deg
	
/*	printf("%lx", tab[tabsize - 1]);	*/
	
	res->t[0] = ((uint64_t) stok->t[0]) & (rho - 1);
	//res->t[0] = tab[0] & (rho - 1);
	counter = 0;
	for(i = 1; i < N; i++)
	{
/*		printf("%016lx", tab[tabsize - 1 - i]);*/
		counter = (counter + 64 - RHO) % RHO;
		res->t[i] = ((((uint64_t) stok->t[i]) << counter) & (rho - 1)) | (((uint64_t) stok->t[i - 1]) >> (64 - counter));
		//res->t[i] = ((tab[i] << counter) & (rho - 1)) | (tab[i - 1] >> (64 - counter));
	}
/*	printf("\n");*/
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);

	free_poly(stok);
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

void __init_tests__(void)
{
	poly a, b, c, c_check, Phisquared;
	init_polys(N, &a, &b, &c, &c_check, &Phisquared, NULL);
	
	const char* A = "9b6afe91a6e17ff3e5b7331bc220a825e6bbe48687ca568a0873800b48471d633375";
	poly converted;
	printf("%p\n", converted);
	init_poly(N, &converted);
	printf("%p deg: %d %p\n", converted, converted->deg, converted->t);
	convert_string_to_amns(converted, A);
	printf("%p deg: %d %p\n", converted, converted->deg, converted->t);
	//printf("%lx\n", converted->t[N - 1]);
	printf("%s\n", A);
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

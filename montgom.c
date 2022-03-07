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

static inline void mp_print(poly P)
{
/*	if(P->t[P->deg - 1] & 0x8000000000000000)*/
/*	{*/
/*		printf("-%lx", -P->t[P->deg - 1]);*/
/*		for(register uint16_t i = 1; i < P->deg; i++)*/
/*			printf("%016lx", -P->t[P->deg - 1 - i]);*/
/*	}*/
/*	else*/
/*	{*/
		printf("%lx", P->t[P->deg - 1]);
		for(register uint16_t i = 1; i < P->deg; i++)
			printf("%016lx", P->t[P->deg - 1 - i]);
	//}
	printf("\n"); 
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

static inline void mp_reduce(restrict poly A)
{
	while(A->deg > 1 && A->t[A->deg - 1] == 0) --A->deg;
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
	
	mp_reduce(*res);
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
	printf("%lx", stok->t[N - 1]);
	for(i = 1; i < N; i++)
		printf("%016lx", stok->t[N - 1 - i]);
	printf("\n");
	if(stok->deg != N)
	{
		printf("ERROR: polynomial degree too high in given number for conversion\n");
		goto end;
	}
	
	
	res->t[0] = ((uint64_t) stok->t[0]) & (rho - 1);
	counter = 0;
	for(i = 1; i < N; i++)
	{
		counter = (counter + 64 - RHO) % RHO;
		res->t[i] = ((((uint64_t) stok->t[i]) << counter) & (rho - 1)) | (((uint64_t) stok->t[i - 1]) >> (64 - counter));
	}
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);

end:
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

static inline void mp_copy(restrict poly* A, restrict const poly B)
{
	// Copy B into A.
	if((*A)->deg < B->deg)
	{
		free_poly(*A);
		init_poly(B->deg, A);
	}
	
	(*A)->deg = B->deg;
	for(register uint16_t i = 0; i < B->deg; i++)
		(*A)->t[i] = B->t[i];
}

static inline int8_t mp_ucomp(restrict poly A, restrict poly B)
{
	// compares A and B. Returns 0 if equal, 1 if A > B and -1 if B < A.
	// doesn't check sign.
	
	register int32_t i;
	
	mp_reduce(A);
	mp_reduce(B);

	if(A->deg > B->deg)
		return 1;
	else
	{
		if(B->deg > A->deg)
			return -1;
		else
		{			
			for(i = A->deg - 1; i > -1; i--)
			{
				if(((uint64_t)A->t[i]) > ((uint64_t)B->t[i]))
					return 1;
				if(((uint64_t)A->t[i]) < ((uint64_t)B->t[i]))
					return -1;
			}
			return 0;			
		}
	}
}

static inline int8_t mp_comp(restrict poly A, restrict poly B)
{
	// compares A and B. Returns 0 if equal, 1 if A > B and -1 if B < A.
	// checks sign.
	
	mp_reduce(A);
	mp_reduce(B);

	if((A->t[A->deg - 1] & 0x8000000000000000) &&
			(!(B->t[B->deg - 1] & 0x8000000000000000)))
		return -1;
	if((!(A->t[A->deg - 1] & 0x8000000000000000)) &&
			(B->t[B->deg - 1] & 0x8000000000000000))
		return 1;

	return mp_ucomp(A, B);
}

static inline void mp_leftshift(restrict poly* A)
{
	// Shifts *A to the left once.
	
	poly aux;
	init_poly(0, &aux);
	
	mp_reduce(*A);
	mp_copy(&aux, *A);
	
	if((*A)->t[(*A)->deg - 1] & 0x8000000000000000)
	{
		free_poly(*A);
		init_poly(aux->deg + 1, A);
		(*A)->t[(*A)->deg - 1] = 1;
	}
	
	for(register uint16_t i = 0; i < aux->deg - 1; i++)
		(*A)->t[aux->deg - 1 - i] = (aux->t[aux->deg - 1 - i] << 1) +
			((aux->t[aux->deg - 2 - i] & 0x8000000000000000) != 0);
	
	(*A)->t[0] = aux->t[0] << 1;
	
	free_poly(aux);
}

static inline void mp_rightshift(restrict poly A)
{
	// Shifts A to the right once.
	
	mp_reduce(A);
	
	for(register uint16_t i = 0; i < A->deg - 1; i++)
		A->t[i] = (((uint64_t)A->t[i]) >> 1) | ((A->t[i + 1] & 1) << 63);
	A->t[A->deg - 1] = (((uint64_t)A->t[A->deg - 1]) >> 1);
	
	mp_reduce(A);
}

static inline void mp_add(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 + op2 in res.

	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg + 1 -
		((op1->t[op1->deg - 1] < 0) || (op2->t[op2->deg - 1] < 0));
	register uint16_t i, j;
	uint64_t stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = ((uint64_t) (*res)->t[i]);
		(*res)->t[i] += op2->t[i];
		
		j = i;
		while(stok > ((uint64_t) (*res)->t[j]) && j < MAXDEG - 1)
		{
			++j;
			stok = ((uint64_t) (*res)->t[j]);
			(*res)->t[j] = ((uint64_t) (*res)->t[j]) + 1;
		}
	}
	
	mp_reduce(*res);
}

static inline void mp_sub(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 - op2 in res.

	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	register uint16_t i, j;
	uint64_t stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	
	for(i = 0; i < op1->deg; i++)
		(*res)->t[i] = op1->t[i];
	
	for(i = 0; i < op2->deg; i++)
	{
		stok = ((uint64_t) (*res)->t[i]);
		(*res)->t[i] -= op2->t[i];
		
		j = i;
		while(stok < ((uint64_t) (*res)->t[j]) && j < MAXDEG - 1)
		{
			++j;
			stok = ((uint64_t) (*res)->t[j]);
			(*res)->t[j] = ((uint64_t) (*res)->t[j]) - 1;
		}
	}
	
	mp_reduce(*res);
}

static inline void mp_mult(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 * op2 into res.

	register uint16_t i, j, k;
	unsigned __int128 R[op1->deg + op2->deg], stok;

	// We check if the degree is high enough. If it isn't we fix the problem.
	if((*res)->deg < op1->deg + op2->deg)
	{
		free_poly(*res);
		init_poly(op1->deg + op2->deg, res);
	}
	
	for(i = 0; i < op1->deg + op2->deg; i++)
		R[i] = 0;

	for(i = 0; i < op1->deg; i++)
	{
		for(j = 0; j < op2->deg; j++)
		{
			k = i + j;
			stok = R[k];
			R[k] += (unsigned __int128) ((uint64_t) op1->t[i]) * ((uint64_t)op2->t[j]);
			while(stok > R[k])
			{
				++k;
				stok = R[k];
				R[k] += 1;
			}
		}
	}
	
	(*res)->t[0] = R[0];
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
		(*res)->t[i] = R[i];
	}
	
	mp_reduce(*res);
}

static inline void mp_mod(restrict poly* res, restrict const poly op1, restrict const poly op2)
{
	// Puts the result of op1 % op2 into res.
	
	const uint16_t MAXDEG = op1->deg < op2->deg ? op2->deg : op1->deg;
	poly X, R;
	init_polys(MAXDEG, &X, &R, NULL);
	
	if((*res)->deg < MAXDEG)
	{
		free_poly(*res);
		init_poly(MAXDEG, res);
	}
	
	mp_copy(&X, op2);
	while(mp_ucomp(op1, X) == 1)
		mp_leftshift(&X);
	mp_print(X);
	while(mp_ucomp(op1, X) == -1)
		mp_rightshift(X);
	mp_print(X);
	
	mp_sub(res, op1, X);
	
	while(mp_ucomp(*res, op2) == 1)
	{
		mp_copy(&R, *res);
		while(mp_ucomp(X, *res) == 1)
		{
			mp_rightshift(X);
		}
		mp_sub(res, R, X);
	}
	
	mp_reduce(*res);
	free_polys(X, R, NULL);
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


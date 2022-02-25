#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>


#define RHO 61
#define N 5
#define LAMBDA 2

typedef int64_t poly[N];
poly M = {9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940}, M1 = {7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353}, MLambda = {18892189295192050, -178689719550530, 732756003058628, -9116351660287002, 38462084564469880}, M1Lambda = {-4355481239625866354, 6277017080637277974, -2433945398614454754, 3202559735019019384, 8710962479251732706};

void set_val(poly P, int64_t val, ...)
{
	va_list args;
	
	va_start(args, val);
	for(int16_t i = 0; i < N; i++)
	{
		P[i] = val;
		val = va_arg(args, int64_t);
	}
	va_end(args);
}

static inline void print(const poly P)
{
	printf("[");
	for(int16_t i = 0; i < N - 1; i++)
		printf("%ld, ", P[i]);
	printf("%ld]\n", P[N - 1]);
}

static inline void __print128(register const __int128 Val)
{
	int64_t hi = Val >> 64;
	uint64_t lo = Val;
	if (hi & 0x8000000000000000) printf("-");
	printf("0x%lx%016lx\n", hi & 0x8000000000000000 ? -hi : hi, lo);
}

static inline void mns_mod_mult_ext_red(__int128* R, const poly A,
	const poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction.
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A[i + j] * B[N - j];
		
		R[i] = R[i] * LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * B[i - j];
	}
}

static inline void m_mns_mod_mult_ext_red(__int128* R, const poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A[i + j] * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * M[i - j];
	}
}

static inline void m1_mns_mod_mult_ext_red(__int128* R, const poly A)
{
	// Same as above but with some pre calculations done in the case of M being
	// the second operand.
	
	register uint16_t i, j;
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A[i + j] * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * M1[i - j];
	}
}

static inline void amns_montg_mult(poly S, const poly A,
	const poly B)
{
	// Function that multiplies A by B using the montgomery approach in an
	// amns. Puts the result in res. Needs M a line of the LLL'd base matrix
	// of the set of polynomials of that amns who have gamma as a root such that
	// gcd of M and E is equal to an odd number. M1 is -((M^-1) mod E) mod phi).
	register uint16_t i, j;
	__int128 R[N] = {0};
	uint64_t V[N], V2[N], T[N], T2[N], Q[N];
	
	for(i = 0; i < N; i++)
	{
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) A[i + j] * B[N - j];
		
		R[i] = R[i] * LAMBDA;
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * B[i - j];
		V[i] = R[i];
		V2[i] = (R[i] >> 64);
	}
	
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) V[i + j] * M1Lambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) V[j] * M1[i - j];
		Q[i] = R[i];
	}
	
	for(i = 0; i < N; i++)
	{
		R[i] = 0;
		for(j = 1; j < N - i; j++)
			R[i] += (__int128) Q[i + j] * MLambda[N - j];
		
		for(j = 0; j < i + 1; j++)
			R[i] += (__int128) Q[j] * M[i - j];
		
		T[i] = R[i];
		T[i] = V[i] + T[i];
		T2[i] = (R[i] >> 64);
		S[i] = V2[i] + T2[i] + (T[i] < V[i]);
	}
/*	mns_mod_mult_ext_red(R, A, B);*/
/*	*/
/*	for(int i = 0; i < N; i++)*/
/*	{*/
/*		V[i] = R[i];*/
/*		S[i] = R[i];*/
/*		V2[i] = (R[i] >> 64);*/
/*		//__print128(R[i]);*/
/*		R[i] = 0;*/
/*	}*/
/*	*/
/*	//print(res);*/
/*	m1_mns_mod_mult_ext_red(R, S);*/
/*	*/
/*	for(int i = 0; i < N; i++)*/
/*	{*/
/*		S[i] = R[i];*/
/*		//__print128(R[i]);*/
/*		R[i] = 0;*/
/*	}*/
/*	*/
/*	m_mns_mod_mult_ext_red(R, S);*/
/*	*/
/*	for(int i = 0; i < N; i++)*/
/*	{*/
/*		T[i] = R[i];*/
/*		T2[i] = (R[i] >> 64);*/
/*		*/
/*		T[i] = V[i] + T[i];*/
/*		S[i] = V2[i] + T2[i] + (T[i] < V[i]);*/
/*		//__print128(R[i]);*/
/*	}*/

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
	for(register uint16_t i = 0; i < N; i++)
		P[i] = __modrho(randomint64());
}

void __init_tests__(void)
{
	
	poly a = {3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128};
	poly b = {1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068};
	//poly b = {1, 0, 0, 0, 0};
	poly c;
	//poly c_check, Phisquared = {0, 0, 0, 512, 0};
	
	//amns_montg_mult(c_check, a, Phisquared);
	//amns_montg_mult(c, c_check, b);

	amns_montg_mult(c, a, b);

	print(a);
	print(b);
	//print(c_check);
	print(c);
	
}

void __proof_of_accuracy(void)
{
	time_t seed;
	poly a, b, c;
	srand((unsigned)time(&seed));
	
	for(register int64_t i = 0; i < 100000; i++)
	{
		randpoly(a);
		randpoly(b);
		amns_montg_mult(c, a, b);
		
		print(a);
		print(b);
		print(c);
	}
	
}

int main(void)
{
	__init_tests__();
	//__proof_of_accuracy();
	return 0;
}


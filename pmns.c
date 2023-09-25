#include <stdio.h>
#include <time.h>
#include <omp.h>

#include "pmns.h"
#include "params.h"
#include "utilitymp.h"

#define NBTHREADZ 8
#define MULTITHREA

const uint8_t PDEGREE = N;

#if PHI > 64

const unsigned __int128 PHIMASK = (((__int128) ((1 << (PHI - 64)) - 1) << 64) |
	((uint64_t)-1L));

void mns_mod_mult_ext_red(__int128* restrict R, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	int64_t matrix[2*N - 1];
	
	for(int i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, B, R) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		R[i] = (__int128) A->t[0] * matrix[N - 1 + i];
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N; j++)
			R[i] += (__int128) A->t[j] * matrix[N - 1 - j + i];
	}
}

void m_or_b_mns_mod_mult_ext_red(__int128* restrict R, int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, matrM, R) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		//printf("0x"); __print128(R[i]); printf("\n");
		R[i] = (R[i] >> (PHI - 64)) + (__int128) A[0] * matrM[N - 1 + i];
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N; j++)
			R[i] += (__int128) A[j] * matrM[N - 1 - j + i];
		
	}
}

void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R, __int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	__int128 somme;
	//__int128 soma;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, matrM1, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		somme = 0;
		//soma = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
		{
			//soma += ( A[j] >> (PHI - 64)) * matrM1[N - 1 - j + i];
			somme += matrM1[N - 1 - j + i] * A[j];
		}
		//printf("0x"); __print128(somme >> (PHI - 64)); printf("\n");
		//printf("0x"); __print128(soma); printf("\n\n");
		R[i] = somme >> (PHI - 64);
	}
}

inline void mns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Montgomery method.
	int64_t T[N];
	register uint16_t i;
	
	m1_or_b1_mns_mod_mult_ext_red(T, R);
	
	m_or_b_mns_mod_mult_ext_red(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] < 0);
}

inline void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N];
	
	mns_mod_mult_ext_red(R, A, B);
	
	mns_montg_int_red(res, R);
}

#else

#ifdef LAMBDA
void new_matvec_toeplitz_par(poly B, poly A, __int128 *res, int size)
{

	int64_t *vec = A->t;
 int i, j;
 int64_t *aux;
__int128 aux2;

	int64_t matrix[2*N - 1];
	
	for(i = 0; i < N-1; i++)
	{
		matrix[i + N - 1] = B->t[i];
		matrix[i] = B->t[1 + i] * LAMBDA;
	}
	matrix[2*N - 2] = B->t[N - 1];

//omp_set_num_threads(8); 
#pragma omp parallel for shared(size,matrix,vec,res) private(i,j,aux,aux2) num_threads(NBTHREADZ)

 for (i = 0; i < size; i++)
 {
     aux2 = 0;
     aux = matrix + 382 - i;
     for (j = 0; j < size; j++)
     {
  //       printf("%d %d %ld %ld : ",i,j,matrix[size-1-i+j],vec[j]);
          aux2 += (__int128)*(aux--)*vec[j];
     }
     res[size-1-i] = aux2;
    // printf("\n");
  }
}

void normmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R
	
	__int128 somme;
	
	//int64_t matrix[2*N - 1];
	
	//for(int i = 0; i < N-1; i++)
	//{
	//	matrix[i + N - 1] = B->t[i];
	//	matrix[i] = B->t[1 + i] * LAMBDA;
	//}
	//matrix[2*N - 2] = B->t[N - 1];
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, B, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		//somme = 0;
		//_PRAGMAGCCUNROLLLOOP_
		//for(int j = 0; j < N; j++)
		//	somme += (__int128) A->t[j] * matrix[N - 1 - j + i];
		//R[i] = somme;
		
		somme = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
			somme += (__int128) A->t[i + j] * B->t[N - j];
		
		somme *= LAMBDA;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
			somme += (__int128) A->t[j] * B->t[i - j];
		R[i] = somme;
	}
}

void toeplitzmns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Toeplitz version.
	
	//register uint16_t i;
	
	int64_t matr[2*N - 1];
	
	//_PRAGMAGCCUNROLLLOOP_
	//for(i = 0; i < N-1; i++)
	//{
	//	matr[i] = B->t[i];
	//	matr[i + N] = B->t[N - 1 - i] * LAMBDA;
	//}
	//matr[N - 1] = B->t[N - 1];
	
	for(int i = 0; i < N-1; i++)
	{
		matr[i + N - 1] = B->t[i];
		matr[i] = B->t[1 + i] * LAMBDA;
	}
	matr[2*N - 2] = B->t[N - 1];
	
	//_PRAGMAGCCUNROLLLOOP_
	//for(i = 0; i < N-1; i++)
	//{
	//	matr[N - 1 - i] = B->t[i];
	//	matr[N + i] = B->t[N - 1 - i] * LAMBDA;
	//}
	//matr[0] = B->t[N - 1];
	
	//toeplitz_vm3x3_192x192(R, A->t, matr);
	//toeplitz_vm(R, A->t, matr, N);
}

void karamns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Karatsuba version.
	
	register uint16_t i, j;
	
	__int128 Lo[N/2] = {0}, Hi[N] = {0}, Mid[N] = {0};
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
		{
			Lo[i] += (__int128) A->t[i - j] * B->t[j];
			Hi[i] += (__int128) A->t[N/2 + i - j] * B->t[N/2 + j];
			Mid[i] += (__int128) (A->t[i - j] + A->t[i - j + N/2]) * (B->t[j] + B->t[j + N/2]);
		}
		_PRAGMAGCCUNROLLLOOP_
		for(j = i + 1; j < N/2; j++)
		{
			Hi[i] -= (__int128) A->t[i + N/2 - j] * B->t[j];
			Hi[i + N/2] += (__int128) A->t[N + i - j] * B->t[N/2 + j];
			Mid[i + N/2] += (__int128) (A->t[i + N/2 - j] + A->t[i - j + N]) * (B->t[j] + B->t[j + N/2]);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
	{
		R[i + N/2] +=  Hi[i + N/2] * LAMBDA + Mid[i] - Lo[i] - Hi[i];
		R[i] += Lo[i] + (Mid[i + N/2] + Hi[i] - Hi[i + N/2]) * LAMBDA;
	}
}

#ifdef M_or_B_is_M
void normm_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	const int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	__int128 somme;
	
	//_PRAGMAGCCUNROLLLOOP_
	//omp_set_num_threads(8);
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, M, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		somme = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			somme += (__int128) A[j] * matrM[N - 1 - j + i];
		//_PRAGMAGCCUNROLLLOOP_
		//for(int j = 1; j < N - i; j++)
		//	somme += (__int128) (A[i + j]) * M[N - j];
		
		//somme *= LAMBDA;
		
		//_PRAGMAGCCUNROLLLOOP_
		//for(int j = 0; j < i + 1; j++)
		//	somme += (__int128) (A[j]) * M[i - j];
		
		R[i] += somme;
	}
}

void toeplitzm_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Toeplitz version.
	
	//Mtoeplitz_vm3x3_192x192(R, A);//, matrM);
	//nonovtoeplitz_vm3x3_192x192(R, A, matrM);
	//M192_toep_3x3(R, A);
}

void karam_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Karatsuba version.
	
	register uint16_t i, j;
	
	__int128 Lo[N/2] = {0}, Hi[N] = {0}, Mid[N] = {0}, tmp;
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
		{
			Lo[i] += (__int128) A[i - j] * M[j];
			Hi[i] += (__int128) A[N/2 + i - j] * M[N/2 + j];
			tmp = (M[j] + M[j + N/2]);
			Mid[i] += (__int128) A[i - j]  * tmp;
			Mid[i] += (__int128) A[i - j + N/2] * tmp;
		}
		_PRAGMAGCCUNROLLLOOP_
		for(j = i + 1; j < N/2; j++)
		{
			Hi[i] -= (__int128) A[i + N/2 - j] * M[j];
			Hi[i + N/2] += (__int128) A[N + i - j] * M[N/2 + j];
			tmp = (M[j] + M[j + N/2]);
			Mid[i + N/2] += (__int128) A[i + N/2 - j] * tmp;
			Mid[i + N/2] += (__int128) A[i - j + N] * tmp;
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
		R[i + N/2] +=  Hi[i + N/2] * LAMBDA + Mid[i] -  Lo[i] -  Hi[i];
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
		R[i] += Lo[i] + (Mid[i + N/2] + Hi[i] - Hi[i + N/2]) * LAMBDA;
}

void normm1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	__int128 somme;
	
	//_PRAGMAGCCUNROLLLOOP_
	//omp_set_num_threads(8);
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, M1, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		R[i] = 0;
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < N; j++)
			R[i] += (int64_t) A[j] * matrM1[N - 1 - j + i];
		//_PRAGMAGCCUNROLLLOOP_
		//for(int j = 1; j < N - i; j++)
		//	somme += ((uint64_t)A[i + j]) * M1[N - j];
		
		//somme *= LAMBDA;
		
		//_PRAGMAGCCUNROLLLOOP_
		//for(int j = 0; j < i + 1; j++)
		//	somme += ((uint64_t)A[j]) * M1[i - j];
		
		//R[i] = somme;
	}
}

void toeplitzm1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Toeplitz version.
	
	//M1toeplitz_vm3x3_192x192(R, A);//, matrM1);
}

void karam1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	// Karatsuba version.
	
	register uint16_t i, j;
	
	uint64_t Lo[N/2] = {0}, Hi[N] = {0}, Mid[N] = {0};
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
	{
		_PRAGMAGCCUNROLLLOOP_
		for(j = 0; j < i + 1; j++)
		{
			Lo[i] += A[i - j] * M1[j];
			Hi[i] += A[N/2 + i - j] * M1[N/2 + j];
			Mid[i] += (A[i - j] + A[i - j + N/2]) * (M1[j] + M1[j + N/2]);
		}
		_PRAGMAGCCUNROLLLOOP_
		for(j = i + 1; j < N/2; j++)
		{
			Hi[i] -= A[i + N/2 - j] * M1[j];
			Hi[i + N/2] += A[N + i - j] * M1[N/2 + j];
			Mid[i + N/2] += (A[i + N/2 - j] + A[i - j + N]) * (M1[j] + M1[j + N/2]);
		}
	}
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
		R[i + N/2] = Hi[i + N/2] * LAMBDA + Mid[i] - Lo[i] - Hi[i];
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N/2; i++)
		R[i] = Lo[i] + (Mid[i + N/2] + Hi[i] -  Hi[i + N/2]) * LAMBDA;
}

void new_matvec_toeplitz_parplus(int64_t *vec, __int128 *res, int size)
{

 int i, j;
 int64_t *aux;
__int128 aux2;
const int64_t* matrix = &(matrM[0]);


//omp_set_num_threads(8); 
#pragma omp parallel for shared(size,matrix,vec,res) private(i,j,aux,aux2) num_threads(NBTHREADZ)

 for (i = 0; i < size; i++)
 {
     aux2 = 0;
     aux = matrix + 382 - i;
     for (j = 0; j < size; j++)
     {
  //       printf("%d %d %ld %ld : ",i,j,matrix[size-1-i+j],vec[j]);
          aux2 += (__int128)*(aux--)*vec[j];
     }
     res[size-1-i] += aux2;
    // printf("\n");
  }
}

void new_matvec_toeplitz_par64(const int64_t *matrix, __int128 *vec, int64_t *res, int size)
{

 int i, j;
 int64_t *aux;
uint64_t aux2;

//omp_set_num_threads(8); 
#pragma omp parallel for shared(size,matrix,vec,res) private(i,j,aux,aux2) num_threads(NBTHREADZ)

 for (i = 0; i < size; i++)
 {
     aux2 = 0;
     aux = matrix + 382 - i;
     for (j = 0; j < size; j++)
     {
  //       printf("%d %d %ld %ld : ",i,j,matrix[size-1-i+j],vec[j]);
          aux2 += (uint64_t)*(aux--)*vec[j];
     }
     res[size-1-i] = aux2;
    // printf("\n");
  }
}

void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	//toeplitzm_or_b_mns_mod_mult_ext_red(R, A);
	//Mtoeplitz_vm3x3_189x189(R, A);
	//Mtoeplitz_vm3x3_84x84(R, A);
	//Mtoeplitz_vm3x3_42x42(R, A);
	//Mtoeplitz_vm_40x40(R, A);
	//Mtoeplitz_vm_32x32(R, A);
	//Mtoeplitz_vm20x20(R, A);
	normm_or_b_mns_mod_mult_ext_red(R, A);
	//new_matvec_toeplitz_parplus(A, R, N);
}

void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	//toeplitzm1_or_b1_mns_mod_mult_ext_red(R, A);
	//M1toeplitz_vm3x3_189x189(R, A);
	//M1toeplitz_vm3x3_84x84(R, A);
	//M1toeplitz_vm3x3_42x42(R, A);
	//M1toeplitz_vm_40x40(R, A);
	//M1toeplitz_vm_32x32(R, A);
	//M1toeplitz_vm20x20(R, A);
	normm1_or_b1_mns_mod_mult_ext_red(R, A);
	//new_matvec_toeplitz_par64(matrM1, A, R, N);
}

#endif

#endif

inline void mns_mod_mult_ext_red(__int128* restrict R, const restrict poly A,
	const restrict poly B)
{
	normmns_mod_mult_ext_red(R, A, B);
	//toeplitzmns_mod_mult_ext_red(R, A, B);
	//int64_t matr[2*N - 1];
	//for(int i = 0; i < N-1; i++) 
	//{
	//	matr[i + N - 1] = B->t[i];
	//	matr[i] = B->t[1 + i] * LAMBDA;
	//}
	//matr[2*N - 2] = B->t[N - 1];
	//toeplitz_vm3x3_189x189(R, A->t, matr);
	//toeplitz_vm3x3_84x84(R, A->t, matr);
	//toeplitz_vm3x3_42x42(R, A->t, matr);
	//toeplitz_vm_40x40(R, A->t, matr);
	//toeplitz_vm_32x32(R, A->t, matr);
	//toeplitz_vm_20x20(R, A->t, matr);
	//new_matvec_toeplitz_par(B, A, R, N);
}

#ifdef LENEXTPOLY

inline void mns_mod_mult_ext_red(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
	// Function that multiplies A by B and applies external reduction using
	// E(X) an irreducible polynomial used for reduction. Result in R
	__int128 T;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, B, R) private(T) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
			T += (__int128) A->t[i + j] * B->t[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
			R[i] += (__int128) A->t[j] * B->t[i - j];
	}
}

#ifdef M_or_B_is_M

void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Function that multiplies A by M and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	__int128 T;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, M, R) private(T) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
			T += (__int128) A[i + j] * M[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
			R[i] += (__int128) A[j] * M[i - j];
	}
}

void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Function that multiplies A by M1 and applies external reduction using
	// E(X) = X^n - lambda as a polynomial used for reduction. Result in R.
	
	__int128 T;
	
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, M1, R) private(T) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		T = 0;
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 1; j < N - i; j++)
			T += ((uint64_t)A[i + j]) * M1[N - j];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int k = 0; (k < LENEXTPOLY) && (k < N-1); k++)
			R[i + k] += T * EXTPOLY[k];
		
		_PRAGMAGCCUNROLLLOOP_
		for(int j = 0; j < i + 1; j++)
			R[i] += ((uint64_t)A[j]) * M1[i - j];
	}
}

#endif

#endif


#ifdef M_or_B_is_B

void m_or_b_mns_mod_mult_ext_red(__int128* restrict R,
	int64_t* restrict A)
{
	// Vector-Matrix multiplication between A and B, result in R.
	
	__int128 somme;
	
	//_PRAGMAGCCUNROLLLOOP_
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, B, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		somme = (__int128)B[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (__int128)B[j][i]*A[j];
		
		R[i] += somme;
	}
}

void m1_or_b1_mns_mod_mult_ext_red(int64_t* restrict R,
	__int128* restrict A)
{
	// Vector-Matrix multiplication between A and B1, result in R.
	
	uint64_t somme;
	
	//_PRAGMAGCCUNROLLLOOP_
	#ifdef MULTITHREAD
	#pragma omp parallel for shared(A, B1, R) private(somme) num_threads(NBTHREADZ)
	#else
	_PRAGMAGCCUNROLLLOOP_
	#endif
	for(int i = 0; i < N; i++)
	{
		somme = (uint64_t)B1[0][i]*A[0];
		for (int j = 1; j < N; j++)
			somme += (uint64_t)B1[j][i]*A[j];
		R[i] = somme;
	}
}

#endif

inline void mns_montg_int_red(restrict poly res, __int128* restrict R)
{
	// Internal reduction of R via the Montgomery method.
	int64_t T[N];
	register uint16_t i;
	
	m1_or_b1_mns_mod_mult_ext_red(T, R);
	
	m_or_b_mns_mod_mult_ext_red(R, T);
	
	_PRAGMAGCCUNROLLLOOP_
	for(i = 0; i < N; i++)
		res->t[i] = (R[i] >> 64) + ((int64_t) R[i] != 0);
}

inline void amns_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B)
{
	// Function that multiplies A by B using the Montgomery approach in an
	// amns. Puts the result in res. A and B have to be in the system and res
	// will be in the pmns also such that if A(gamma) = a * phi mod p and 
	// B(gamma) = b * phi mod p then res(gamma) = a * b * phi mod p
	
	__int128 R[N] = {0};
	
	mns_mod_mult_ext_red(R, A, B);

	mns_montg_int_red(res, R);
}

#endif



void amns_rtl_sqandmult(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Right to Left version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp, num, stok;
	init_polys(N, &tmp, &num, &stok, NULL);
	
	poly_copy(res, &__theta__);
	poly_copy(stok, base);
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[i];
		// We do it in pairs to avoid the use of poly_copy
		for(j = 0; j < 32; j++)
			{
				if(aux & 1)
					amns_montg_mult(tmp, res, stok);
				else
					amns_montg_mult(tmp, res, &__theta__);
				amns_montg_mult(num, stok, stok);
				aux >>= 1;
				if(aux & 1)
					amns_montg_mult(res, tmp, num);
				else
					amns_montg_mult(res, tmp, &__theta__);
				amns_montg_mult(stok, num, num);
				aux >>= 1;
			}
	}
	// If our exponent has an odd number of bits we do this one time too many.
	// For this reason unless we find another way, Left to Right is preferred.
	aux = exponent->t[exponent->deg - 1];
	while(aux)
	{
		if(aux & 1)
			amns_montg_mult(tmp, res, stok);
		else
			amns_montg_mult(tmp, res, &__theta__);
		amns_montg_mult(num, stok, stok);
		aux >>= 1;
		if(aux & 1)
			amns_montg_mult(res, tmp, num);
		else
			amns_montg_mult(res, tmp, &__theta__);
		amns_montg_mult(stok, num, num);
		aux >>= 1;
	}
	
	free_polys(tmp, num, stok, NULL);
}

void amns_ltr_sqandmult(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Left to Right version.
	
	register uint16_t i;
	register uint8_t j;
	register uint64_t aux;
	
	poly tmp;
	init_poly(N, &tmp);
	
	poly_copy(res, &__theta__);
	
	aux = exponent->t[exponent->deg - 1];
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		amns_montg_mult(tmp, res, res);
		if(aux & (1ULL << (63 - j)))
			amns_montg_mult(res, tmp, base);
		else
			amns_montg_mult(res, tmp, &__theta__);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			amns_montg_mult(tmp, res, res);
			if(aux & (1ULL << (63 - j)))
				amns_montg_mult(res, tmp, base);
			else
				amns_montg_mult(res, tmp, &__theta__);
		}
	}
	
	free_poly(tmp);
}

void amns_montg_ladder(restrict poly res, const restrict poly base,
	const restrict mpnum exponent)
{
	// Function for fast exponentiation using the Montgomery ladder.
	// Returns base^exponent % p.
	
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
	for(j = __builtin_clz(aux); j < 64; j++)
	{
		b = (aux & (1ULL << (63 - j))) >> (63 - j);
		amns_montg_mult(R[1 - b], R[1 - b], R[b]);
		amns_montg_mult(R[b], R[b], R[b]);
	}
	
	for(i = 0; i < exponent->deg - 1; i++)
	{
		aux = exponent->t[exponent->deg - 2 - i];
		for(j = 0; j < 64; j++)
		{
			b = (aux & (1ULL << (63 - j))) >> (63 - j);
			amns_montg_mult(R[1 - b], R[1 - b], R[b]);
			amns_montg_mult(R[b], R[b], R[b]);
		}
	}
	
	free_poly(tmp);
}

inline void amns_sqandmult(restrict poly res, const char* restrict base,
	const char* restrict exponent)
{
	// Function for fast exponentiation using the square and multiply algorithm.
	// Returns base^exponent % p. Uses Left to Right square and multiply.
	
	poly pbase;
	mpnum mpexponent;
	init_poly(N, &pbase);
	init_mpnum(1, &mpexponent);
	
	convert_string_to_amns(pbase, base);
	convert_string_to_multipre(&mpexponent, exponent);
	
	amns_ltr_sqandmult(res, pbase, mpexponent);
	
	free_poly(pbase);
	free_mpnum(mpexponent);
}

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

int64_t __modrho(int64_t param)
{
	return param & ((1ULL<<RHO) - 1);
}

void randpoly(poly P)
{
	// Function that generates a random polynomial within our PMNS system.
	
	for(register uint16_t i = 0; i < P->deg; i++)
		P->t[i] = __modrho(randomint64()) * (1 + (rand() & 1) * -2);
}

void convert_string_to_amns(restrict poly res, const char* string)
{
	// Function that converts a hexadecimal number given as a string into a
	// polynomial in our representation system (we multiply it by PHI in the
	// process).
	
	uint8_t counter;
	register uint16_t i, j;
	const uint64_t rho = (1ULL<<RHO);
	__int128 R[N] = {0};
	mpnum stok, tmp;
	init_mpnums(N, &stok, &tmp, NULL);
	
	convert_string_to_multipre(&tmp, string);
	
	mp_mod(&stok, tmp, &__P__);
	
	for(i = 1; i < N; i++)
		res->t[i] = 0;
	
	counter = 0;
	res->t[0] = stok->t[0] & (rho - 1);
	for(i = 1; i < stok->deg; i++)
	{
		counter = (counter + RHO) % 64; 
		res->t[i] = (stok->t[i-1] >> counter);
		res->t[i] |= (stok->t[i] << (64 - counter));
		res->t[i] &= (rho - 1);
	}
	counter = (counter + RHO) % 64;
	res->t[i] = (stok->t[i-1] >> counter);
	
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			R[j] += (__int128) res->t[i] * __Pi__[i][j];
	
	mns_montg_int_red(res, R);
	
	free_mpnums(stok, tmp, NULL);
}

void convert_amns_to_multipre(restrict mpnum* res, const restrict poly P)
{
	// Function that converts out of the AMNS system and into a multiprecision
	// number.
	
	register uint16_t i;
	poly a;
	mpnum aux, ag, tmp, tmp2;
	__int128 Quite[N];
	
	init_poly(N, &a);
	init_mpnums(N, &ag, &tmp, &tmp2, NULL);
	init_mpnum(1, &aux);
	if((*res)->deg < N)
	{
		free_mpnum(*res);
		init_mpnum(N, res);
	}
	(*res)->deg = N;
	for(i = 1; i < N; i++)
	{
		(*res)->t[i] = 0;
		Quite[i] = (__int128) P->t[i];
	}
	Quite[0] = (__int128) P->t[0];
	
	mns_montg_int_red(a, Quite);
	
	// We use Horner to get out of the PMNS for now as the other method has issues.
	(*res)->sign = 1 - 2 * (a->t[N - 1] < 0);
	(*res)->t[0] = a->t[N - 1] * (*res)->sign;
	for(i = 0; i < N - 1; i++)
	{
		mp_mult(&tmp, *res, &Gi[0]);
		aux->sign = 1 - 2 * (a->t[N - 2 - i] < 0);
		aux->t[0] = a->t[N - 2 - i] * aux->sign;
		mp_add(&tmp2, tmp, aux);
		mp_mod(res, tmp2, &__P__);
	}
	
	//(*res)->t[0] = a->t[0];
	//for(i = 1; i < N; i++)
	//{
	//	aux->sign = 1 - 2 * (a->t[i] < 0);
	//	aux->t[0] = a->t[i] * aux->sign;
	//	mp_mult(&ag, aux, &Gi[i - 1]);
	//	mp_copy(&tmp, *res);
	//	mp_add(res, tmp, ag);
	//}
	
	//mp_copy(&tmp, *res);
	//mp_mod(res, tmp, &__P__);
	
	free_mpnums(aux, ag, tmp, tmp2, NULL);
	free_poly(a);
}

void __multchecks__(char* nbmults)
{
	// Used as a debug tool to see if the PMNS correctly gives us the proper
	// results with a few random values.
	poly a, b, c;
	init_polys(N, &a, &b, &c, NULL);
	int64_t seed;
	
	register uint64_t cap = 100;
	
	if(nbmults[0] != '\0')
		cap = atoll(nbmults);
	
	srand((unsigned) (time(&seed)));
	
	for(register uint64_t i = 0; i < cap; i++)
	{
		randpoly(a);
		randpoly(b);
		poly_print(a);
		poly_print(b);
		amns_montg_mult(c, a, b);
		poly_print(c);
	}
	
	free_polys(a, b, c, NULL);
}

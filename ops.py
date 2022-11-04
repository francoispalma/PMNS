"""
Module for various miscellaneous needed operations for either generation
purposes or checking/prototyping the C code.
"""

from random import randrange
from sage.all import matrix
#from pyparams128 import p, n, gamma, lam, rho, M_or_B, M1_or_B1, phi
#import convert

def horner_modulo(Poly: list, X: int, modulo: int) -> int:
	"""Polynomial evaluation function using Horner's algorithm, applying a
	modulo.
	
	PARAMETERS
	Poly: The polynomial in list form with the 0th degree in 0th index
	X: The value to evaluate the polynomial in
	modulo: The value to apply a modulo for during the process
	
	For our purposes we use this to convert out of the PMNS by evaluating in
	gamma modulo p.
	"""
	sum_ = 0
	for i in range(len(Poly) - 1, -1, -1):
		sum_ = sum_ * X
		sum_ = sum_ + Poly[i]
	return int(sum_ % modulo)

# Probably removing this soon.
def amns_mod_mult(A, B, n, lam):
	if type(lam) == int:
		return pmns_mod_mult(A, B, n, [lam])
	else:
		return pmns_mod_mult(A, B, n, lam)
#	R = [0] * n
#	for i in range(n):
#		for j in range(1, n - i):
#			R[i] += A[i + j] * B[n - j]
#		R[i] *= lam
#		for j in range(i + 1):
#			R[i] += A[j] * B[i - j]
#	return R

# Old version using polynomials, probably getting removed later.
def montgomery_like_coefficient_reduction(V, n, lam, phi, M, M1):
	Q = [int(int(V[i]) & (phi - 1)) for i in range(n)]
	Q = amns_mod_mult(Q, M1, n, lam)
	Q = [int(int(Q[i]) & (phi - 1)) for i in range(n)]
	T = amns_mod_mult(Q, M, n, lam)
	S = [int(int(int(V[i]) + int(T[i])) >> (phi.bit_length() - 1)) for i in range(n)]
	return S

# Old version using polynomials, probably getting removed later.
def amns_montg_mult(A, B, n, lam, phi, M, M1):
	return montgomery_like_coefficient_reduction(amns_mod_mult(A, B, n, lam),
		n, lam, phi, M, M1)

def montgomery_like_coefficient_reduction_base(A: list, phi: int, B: list,
		B1: list) -> list:
	"""Function to apply an internal reduction to A using Montgomery's method.
	This version uses the basis matrix instead of polynomials.
	
	PARAMETERS
	A: The polynomial in list form we are applying a reduction to
	phi: The size of our memory registers, typically 2**64 or 2**128
	B: The basis matrix of polynomials P such that P(gamma) = 0 mod p
	B1: This matrix is equal to the opposite of the inverse of B mod phi
	
	Returns a polynomial in list form such that its evaluation in gamma is
	equal to A's divided by phi. We guarantee that each coefficient is
	smaller by a factor of phi compared to A.
	"""
	mA = matrix(A)
	mB = matrix(B)
	mB1 = matrix(B1)
	nmA = (((mA % phi) * mB1) % phi) * mB
	res = list(((mA - nmA) / phi)[0])
	return res

def amns_montg_mult_base(vA, vB, n, lam, phi, B, B1):
	return montgomery_like_coefficient_reduction_base(amns_mod_mult(vA, vB, n, lam),
		phi, B, B1)

def pmns_mod_mult(A, B, n, E):
	R = [0] * n
	for i in range(n):
		T = 0
		for j in range(1, n - i):
			T += A[i + j] * B[n - j]
		for k in range(len(E)):
			if k >= n - i:
				break
			R[i + k] += T * E[k]
		for j in range(i + 1):
			R[i] += A[j] * B[i - j]
	return R

def pmns_montg_mult_base(vA, vB, n, lam, phi, B, B1):
	return montgomery_like_coefficient_reduction_base(pmns_mod_mult(vA, vB, n, lam),
		phi, B, B1)

#if __name__ == "__main__":
#	phin = pow(phi, n, p)
#	a = randrange(int(p**0.5), p)
#	b = randrange(int(p**0.5), p)
#	c = a * b % p
#	aargs = (p, n, lam, phi, M_or_B, M1_or_B1)
#	grp = (phi, M_or_B, M1_or_B1)
#	A = convert.montgomery_convert_to_mns_base(a, *aargs, phin)
#	B = convert.montgomery_convert_to_mns_base(b, *aargs, phin)
#	C = montgomery_like_coefficient_reduction_base(pmns_montg_mult_base(A, B, n, lam, *grp), *grp)
#	cc = horner_modulo(C, gamma, p)
#	print(c == cc)

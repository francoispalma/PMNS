from random import randrange
from time import process_time
from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice

from ops import montgomery_like_coefficient_reduction, montgomery_like_coefficient_reduction_base, horner_modulo, naive_convert_to_mns, amns_montg_mult, list_to_poly
from pyparams128 import p, n, gamma, rho, lam, M_or_B as M, M1_or_B1 as M1, phi

def norme_infinie(L):
	return max([abs(elem) for elem in L])

def naive_convert_to_mns_and_red(a, p, n, gamma, rho, lam, phi, M, M1, phisquared):
	A = naive_convert_to_mns(a, p, n, gamma, rho, lam)
	U = amns_montg_mult(A, phisquared, p, n, gamma, rho, lam, phi, M, M1)
	return U

def montgomery_convert_to_mns(a, p, n, gamma, rho, lam, phi, M, M1, tau):
	alpha = (a * tau) % p
	A = [0] * n
	A[0] = alpha
	for _ in range(n - 1):
		A = montgomery_like_coefficient_reduction(A, p, n, gamma, rho, lam, phi, M, M1)
	return A

def rho_div_convert_to_mns(a, p, n, gamma, rho, lam, phi, M, M1, Pi):
	t = convert_to_rho_base(a, n, rho)
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
	A = montgomery_like_coefficient_reduction(U, p, n, gamma, rho, lam, phi, M, M1)
	return A

def convert_to_rho_base(a, n, rho):
	t = [0] * n
	a1 = a
	RHO = rho.bit_length() - 1
	for i in range(n):
		t[i] = a1 & (rho - 1)
		a1 >>= RHO
	return t

def montgomery_convert_to_mns_base(a, p, n, gamma, rho, lam, phi, B, B1, tau):
	alpha = (a * tau) % p
	A = [0] * n
	A[0] = alpha
	for _ in range(n - 1):
		A = montgomery_like_coefficient_reduction_base(A, p, n, gamma, rho, lam, phi, B, B1)
	return A

def rho_div_convert_to_mns_base(a, p, n, gamma, rho, lam, phi, B, B1, Pi):
	t = convert_to_rho_base(a, n, rho)
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
	A = montgomery_like_coefficient_reduction_base(U, p, n, gamma, rho, lam, phi, B, B1)
	return A


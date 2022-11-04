from random import randrange
from time import process_time
from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice

from ops import montgomery_like_coefficient_reduction, montgomery_like_coefficient_reduction_base, horner_modulo, amns_montg_mult
from pyparams128 import p, n, gamma, rho, lam, M_or_B as M, M1_or_B1 as M1, phi

def norme_infinie(L):
	return max([abs(elem) for elem in L])

def montgomery_convert_to_mns(a, p, n, lam, phi, M, M1, tau):
	alpha = (a * tau) % p
	A = [0] * n
	A[0] = alpha
	for _ in range(n - 1):
		A = montgomery_like_coefficient_reduction(A, n, lam, phi, M, M1)
	return A

def rho_div_convert_to_mns(a, n, rho, lam, phi, M, M1, Pi):
	t = convert_to_rho_base(a, n, rho)
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
	A = montgomery_like_coefficient_reduction(U, n, lam, phi, M, M1)
	return A

def convert_to_rho_base(a, n, rho):
	t = [0] * n
	a1 = a
	RHO = rho.bit_length() - 1
	for i in range(n):
		t[i] = a1 & (rho - 1)
		a1 >>= RHO
	return t

def montgomery_convert_to_mns_base(a, p, n, lam, phi, B, B1, tau):
	alpha = (a * tau) % p
	A = [0] * n
	A[0] = alpha
	for _ in range(n - 1):
		A = montgomery_like_coefficient_reduction_base(A, phi, B, B1)
	return A

def rho_div_convert_to_mns_base(a, n, rho, lam, phi, B, B1, Pi):
	t = convert_to_rho_base(a, n, rho)
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
	A = montgomery_like_coefficient_reduction_base(U, phi, B, B1)
	return A

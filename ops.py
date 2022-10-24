from random import randrange
from time import process_time
from sage.all import next_prime, factor, matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from copy import deepcopy

def external_reduction(C, n, lam):
	R = [0] * n
	for i in range(n - 1):
		R[i] = C[i] + lam * C[n + i]
	R[n - 1] = C[n - 1]
	return R

def internal_reduction(V, n, k, p, gamma, M):
	U = [0] * n
	W = [0] * n
	for i in range(n):
		U[i] = V[i] >> k
		W[i] = V[i] & ((1<<k) - 1)
		assert (U[i] << k) + W[i] == V[i]
	for i in range(n):
		for j in range(n):
			W[i] = (W[i] + U[j] * M[j][i]) % p
		W[i] %= p
	assert horner_modulo(V, gamma, p) == horner_modulo(W, gamma, p)
	return W

def naive_convert_to_mns(val, p, n, gamma, rho, lam=0):
	gamms = [int(pow(gamma, i, p)) for i in range(n)]
	mgamms = [(p - elem) % p for elem in gamms]
	R = [0] * n
	left = [i for i in range(n - 1, -1, -1)]
	while left:
		orlen = len(left)
		for elem in left:
			if val // gamms[elem] != 0 and -rho < val // gamms[elem] < rho:
				R[elem] = val // gamms[elem]
				val = (val - R[elem] * gamms[elem]) % p
				left.remove(elem)
				break
			elif val // mgamms[elem] != 0 and -rho < val // mgamms[elem] < rho:
				R[elem] = - (val // mgamms[elem])
				val = (val - R[elem] * gamms[elem]) % p
				left.remove(elem)
				break
		if len(left) == orlen:
			print("Error:", val)
			break
	return R

def horner_modulo(Poly, X, modulo):
	sum_ = 0
	for i in range(len(Poly) - 1, -1, -1):
		sum_ = sum_ * X
		sum_ = sum_ + Poly[i]
	return int(sum_ % modulo)

def list_to_poly(L):
	rstr = str(L[0])
	for i in range(1, len(L)):
		rstr += " + (" + str(L[i]) + ") * X" + ("^" + str(i) if i != 1 else "")
	return rstr

def mns_mod_mult(A, B, p, n, gamma, rho, lam):
	R = [0] * n
	for i in range(n):
		for j in range(1, n - i):
			R[i] += A[i + j] * B[n - j]
		R[i] *= lam
		for j in range(i + 1):
			R[i] += A[j] * B[i - j]
	return R

def montgomery_like_coefficient_reduction(V, p, n, gamma, rho, lam, phi, M, M1):
	Q = [int(int(V[i]) & (phi - 1)) for i in range(n)]
	Q = mns_mod_mult(Q, M1, p, n, gamma, rho, lam)
	Q = [int(int(Q[i]) & (phi - 1)) for i in range(n)]
	T = mns_mod_mult(Q, M, p, n, gamma, rho, lam)
	S = [int(int(int(V[i]) + int(T[i])) >> (phi.bit_length() - 1)) for i in range(n)]
	return S

def amns_montg_mult(A, B, p, n, gamma, rho, lam, phi, M, M1):
	amns = (p, n, gamma, rho, lam)
	return montgomery_like_coefficient_reduction(mns_mod_mult(A, B, *amns),
		*amns, phi, M, M1)

def montgomery_like_coefficient_reduction_base(A, p, n, gamma, rho, lam, phi, B, B1):
	mA = matrix(A)
	mB = matrix(B)
	mB1 = matrix(B1)
	nmA = (((mA % phi) * mB1) % phi) * mB
	res = list(((mA - nmA) / phi)[0])
	return res

def amns_montg_mult_base(vA, vB, p, n, gamma, rho, lam, phi, B, B1):
	amns = (p, n, gamma, rho, lam)
	return montgomery_like_coefficient_reduction_base(mns_mod_mult(vA, vB, *amns),
		*amns, phi, B, B1)

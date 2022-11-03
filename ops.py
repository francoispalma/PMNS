from random import randrange
from time import process_time
from sage.all import next_prime, factor, matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from copy import deepcopy
from pyparams import p, n, gamma, lam, rho, M_or_B, M1_or_B1, phi
#import convert

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

def naive_convert_to_mns(val, p, n, gamma, rho, lam=None):
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

def amns_mod_mult(A, B, n, lam):
	R = [0] * n
	for i in range(n):
		for j in range(1, n - i):
			R[i] += A[i + j] * B[n - j]
		R[i] *= lam
		for j in range(i + 1):
			R[i] += A[j] * B[i - j]
	return R

def montgomery_like_coefficient_reduction(V, n, lam, phi, M, M1):
	Q = [int(int(V[i]) & (phi - 1)) for i in range(n)]
	Q = amns_mod_mult(Q, M1, n, lam)
	Q = [int(int(Q[i]) & (phi - 1)) for i in range(n)]
	T = amns_mod_mult(Q, M, n, lam)
	S = [int(int(int(V[i]) + int(T[i])) >> (phi.bit_length() - 1)) for i in range(n)]
	return S

def amns_montg_mult(A, B, n, lam, phi, M, M1):
	return montgomery_like_coefficient_reduction(amns_mod_mult(A, B, n, lam),
		n, lam, phi, M, M1)

def montgomery_like_coefficient_reduction_base(A, phi, B, B1):
	mA = matrix(A)
	mB = matrix(B)
	mB1 = matrix(B1)
	nmA = (((mA % phi) * mB1) % phi) * mB
	res = list(((mA - nmA) / phi)[0])
	return res

def amns_montg_mult_base(vA, vB, n, lam, phi, B, B1):
	return montgomery_like_coefficient_reduction_base(amns_mod_mult(vA, vB, n, lam),
		phi, B, B1)

def pmns_mod_mult(vA, vB, n, E):
	R = [0] * n
	for i in range(n):
		for j in range(1, n - i):
			T = A[i + j] * B[n - j]
			for k in range(len(E)):
				R[i + k] += T * E[k]
		for j in range(i + 1):
			R[i] += A[j] * B[i - j]
	return R

def pmns_montg_mult_base(vA, vB, n, lam, phi, B, B1):
	return montgomery_like_coefficient_reduction_base(pmns_mod_mult(vA, vB, n, lam),
		phi, B, B1)

#if __name__ == "__main__":
#	phin = pow(phi, n, p)
#	a = randrange(p**0.5, p)
#	b = randrange(p**0.5, p)
#	c = a * b % p
#	pmns = (p, n, gamma, rho, lam, phi, M_or_B, M1_or_B1)
#	grp = (phi, M_or_B, M1_or_B1)
#	A = convert.montgomery_convert_to_mns_base(a, *pmns, phin)
#	B = convert.montgomery_convert_to_mns_base(b, *pmns, phin)
#	C = montgomery_like_coefficient_reduction_base(pmns_montg_mult_base(A, B, n, lam, *grp), *grp)
#	cc = horner_modulo(C, gamma, p)
#	print(c == cc)

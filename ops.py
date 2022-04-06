from random import randrange
from time import process_time
from sage.all import next_prime, factor, matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from copy import deepcopy

# nth root
# b = pgcd(p-1, n)
# alpha^((p-1)/b) = 1 (mod p)
# Tonelli-shanks algorithm?
# 5-11th root (6 = 3x2, etc)
# get gamma from lambda since gamma = nth root of lambda X^n = lambda
# because E(X) = X^n - lambda
# gamma root of E(X)
# generate small lambda, gamma < 64 bits
# pgcd(p-1, n) = 1
# => nu + (p-1)v = 1 (Bezout)
# pgcd(p, lambda) = 1
# => lambda ^ (p-1) = 1 mod p (Euler)
# lambda ^ (p-1)v = 1 mod p
# lambda ^ (nu) * lambda ^ (p-1)v = lambda ^ (nu) mod p
# lambda ^ (nu + (p-1)v) = lambda^(nu) mod p
# lambda = lambda^(nu) mod p
# lambda = (lambda^u)^n mod p
# lambda^u nth root of lambda
# gamma = lambda^u mod p
# else pgcd(p, n) = b
# nu + (p-1)v = b (Bezout)
# lambda ^ (nu + (p-1)v) = lambda^(nu) mod p
# lambda ^ b = lambda ^(nu) mod p
# if u = (b-1)k
# lambda^b = lambda^(n(b-1)k)
# lambda = lambda^(n(b-1)k - (b-1))
# lambda = lambda^((b-1)(nk - 1))
# ???

def trivial_modular_add(a, b, p):
	s = a + b
	if s >= p:
		s = s - p
	return s

def omura_modular_add(a, b, c, n):
	s = a + b
	t = s + c
	if t & (1 << n):
		s = t ^ (1 << n)
	return s

def barrett_modular_mult(a, b, p, n, v):
	c = a * b
	u = c >> (n - 1)
	w = u * v
	q = w >> (n + 1)
	s = c - q * p
	if s >= p:
		s = s - p
	if s >= p:
		s = s - p
	return s

def montgomery_modular_mult(a, b, n, mpm1, r):
	c = a * b
	q = c * mpm1 & (r - 1)
	s = (c + q * p) >> n
	if s >= p:
		s = s - p
	return s

def mersenne_modular_mult(a, b, p, n):
	c = a * b
	c1 = c >> n
	c0 = c & ((1 << n) - 1)
	assert c == (c1 << n) + c0
	s = c1 + c0
	if s >= p:
		s = s - p
	return s

def pseudo_mersenne_modular_mult(a, b, p, n, c1):
	c = a * b
	r1 = c >> (3 * n // 2)
	r0 = c & ((1 << (3 * n // 2)) - 1)
	assert c == (r1 << (3 * n // 2)) + r0
	r = r1 * c1 * (1 << (n // 2)) + r0
	s1 = r >> n
	s0 = r & ((1 << n) - 1)
	assert r == (s1 << n) + s0
	s = s1 * c1 + s0
	if s >= p:
		s = s - p
	return s

def left_to_right_square_and_multiply(a, h, p, n):
	s = 1
	for i in range(n - 1, -1, -1):
		s = (s * s) % p
		if h & (1 << i):
			s = (s * a) % p
	return s

def left_to_right_square_and_multiply_always(a, h, p, n):
	R0 = 1
	for i in range(n - 1, -1, -1):
		R0 = (R0 * R0) % p
		R1 = (R0 * a) % p
		R0 = R1 if (h & (1 << i)) else R0
	return R0

def montgomery_ladder(a, h, p, n):
	R = [1, a]
	for i in range(n - 1, -1, -1):
		b = ((1 << i) & h) != 0
		R[1 - b] = (R[0] * R[1]) % p
		R[b] = (R[b] * R[b]) % p
	return R[0]

def extended_euclidan(a, b):
	# au + bv = d
	u1, u2, u3 = 1, 0, a
	v1, v2, v3 = 0, 1, b
	while v3 != 0:
		q = u3//v3
		t1, t2, t3 = v1, v2, v3
		v1, v2, v3 = [a[0] - q * a[1] for a in zip((u1, u2, u3), (v1, v2, v3))]
		u1, u2, u3 = t1, t2, t3
	u, v, d = u1, u2, u3
	return u, v, d

def swaplines(M, n, lig1, lig2):
	for i in range(n):
		M[lig1][i], M[lig2][i] = M[lig2][i], M[lig1][i]

def gauss_pivot(Input):
	"""Returns determinant of input sing gauss pivot.
	"""
	M = deepcopy(Input)
	n = len(M)
	for col in range(n - 1):
		i = col
		while M[i][col] == 0:
			i += 1
		if i != col:
			swaplines(M, n, i, col)
		for i in range(col + 1, n):
			if M[i][col] != 0:
				store = M[i][col]
				for j in range(col, n):
					M[i][j] = M[i][j] * M[col][col] - store * M[col][j]
	product = 1
	for i in range(n):
		product *= M[i][i]
	return product

def get_prime_factors(val):
	L = list(factor(val))
	return [L[i][0] for i in range(len(L))]

def random_mns(l, n, k, maxlampower):
	assert maxlampower >= 0, "Unacceptable power of two"
	while True:
		lam = 1 << randrange
		ksi = [0] * n
		for _ in range(randrange(1, 3)):
			ksi[randrange(n)] = randrange(-1, 2, 2)
		M = [[ksi[(i - j) % (n)] * ((i - j < 0) * lam + (i - j >= 0)) for i in range(n)] for j in range(n)]
		M = [[-M[i][j] if i != j else (1<<k) -M[i][i] for i in range(n)] for j in range(n)]
		det = matrix(M).determinant()
		L = get_prime_factors(det)
		if L[-1] > (1<<(l-1)):
			break
	p = L[-1]
	rho = (1<<(k + 1))
	u, soak, pgcd = extended_euclidan(n, p-1)
	if pgcd == 1:
		gamma = pow(lam, u, p)
	else:
		# TODO: do something here
		gamma = None
	return (p, n, gamma, rho, lam)

def create_mns(l, n, k, lam, ksi):
	M = [[ksi[(i - j) % (n)] * ((i - j < 0) * lam + (i - j >= 0)) for i in range(n)] for j in range(n)]
	Mprime = [[-M[i][j] if i != j else (1<<k) -M[i][i] for i in range(n)] for j in range(n)]
	det = gauss_pivot(Mprime)
	L = get_prime_factors(det)
	if L[-1] < (1<<(l-1)):
		raise ValueError("Bad parameters, no prime factor big enough found.")
	p = L[-1]
	rho = (1<<(k + 1))
	u, soak, pgcd = extended_euclidan(n, p-1)
	if pgcd == 1:
		gamma = pow(lam, u, p)
	else:
		# TODO: do something here
		gamma = None
	return (p, gamma, rho, M)

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

#def internal_convert_to_mns(val, n, k, p, rho, M):
#	V = [0] * n
#	V[0] = val
#	while True:
#		print(V)
#		V = internal_reduction(V, n, k, p, M)
#		flag = True
#		for elem in V:
#			if elem >= rho:
#				flag = False
#				break
#		if flag == True:
#			break
#	return V

def horner_modulo(Poly, X, modulo):
	sum_ = 0
	for i in range(len(Poly) - 1, -1, -1):
		sum_ = sum_ * X
		sum_ = sum_ + Poly[i]
	return int(sum_ % modulo)

def babai_coefficient_reduction(V, p, n, gamma, rho, B, Betoile):
	S = V.copy()
	for i in range(n):
		c = round(
		int((sum(
		[S[j] * Betoile[n - 1 - i][j] for j in range(n)]
		)) % p ) / int(Betoile[n - 1 - i][n - 1 - i])
		)
		S = [(int(S[j]) - c * int(B[n - 1 - i][j])) for j in range(n)]
	S = [S[i] if abs(S[i]) < p else S[i] % p for i in range(n)]
	return S

def barett_like_coefficient_reduction():
	#TODO: do it. Find a method to find a proper M and how to get M^-1 mod E
	# or don't
	pass

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

if __name__ == "__main__":
	n = 512
	sum1 = 0
	sum2 = 0
	sum3 = 0
	for _ in range(1000):
		p = randrange(1 << (n - 1), 1 << n)
		a = randrange(p)
		b = randrange(p)
		c = (1 << n) - p
		c1 = process_time()
		res1 = trivial_modular_add(a, b, p)
		sum1 += process_time() - c1
		c1 = process_time()
		res2 = omura_modular_add(a, b, c, n)
		sum2 += process_time() - c1
		c1 = process_time()
		res3 = (a + b) % p
		sum3 += process_time() - c1
		if res1 != res2 != res3:
			print(res1)
			print(res2)
			print(res3)
			print()
	print("res:\t", sum1, "\t", sum2, "\t", sum3)

	print("\nMult:")
	sum1 = 0
	sum2 = 0
	sum3 = 0
	for _ in range(1000):
		p = randrange(1 << (n - 1), 1 << n) | 1
		a = randrange(p)
		b = randrange(p)
		v = (1 << n) // p
		r = 1 << n
		mpm1 = pow(-p, -1, r)
		c1 = process_time()
		res1 = barrett_modular_mult(a, b, p, n, v)
		sum1 += process_time() - c1
		ar = (a * r) % p
		br = (b * r) % p
		c1 = process_time()
		res2 = montgomery_modular_mult(ar, br, n, mpm1, r)
		res2 = montgomery_modular_mult(res2, 1, n, mpm1, r)
		sum2 += process_time() - c1
		c1 = process_time()
		res3 = (a * b) % p
		sum3 += process_time() - c1
		if res1 != res2 != res3:
			print("res1:", res1)
			print("res2:", res2)
			print("res3:", res3)
			print()
	print("res:\t", sum1, "\t", sum2, "\t", sum3)

	print("\nMersenne mult:")
	sum1 = 0
	sum2 = 0
	for _ in range(1000):
		n = randrange(512, 1024)
		p = (1 << n) - 1
		a = randrange(p)
		b = randrange(p)
		c1 = process_time()
		res1 = mersenne_modular_mult(a, b, p, n)
		sum1 += process_time() - c1
		c1 = process_time()
		res2 = (a * b) % p
		sum2 += process_time() - c1
		if res1 != res2:
			print("res1:", res1)
			print("res2:", res2)
			print()
	print("res:\t", sum1, "\t", sum2)

	print("\nPseudo Mersenne mult:")
	sum1 = 0
	sum2 = 0
	for _ in range(1000):
		n = (randrange(512, 1024))
		c = randrange(1, 1<<(n//2))
		p = (1 << n) - c
		a = randrange(p)
		b = randrange(p)
		c1 = process_time()
		res1 = pseudo_mersenne_modular_mult(a, b, p, n, c)
		sum1 += process_time() - c1
		c1 = process_time()
		res2 = (a * b) % p
		sum2 += process_time() - c1
#		if res1 != res2:
#			print("res1:", res1)
#			print("res2:", res2)
#			print()
	print("res:\t", sum1, "\t", sum2)

	print("\nModular Exp")
	sum1 = 0
	sum2 = 0
	sum3 = 0
	sum4 = 0
	for _ in range(1000):
		n = (randrange(64, 256))
		p = randrange(1 << n)
		a = randrange(p)
		b = randrange(p)
		c1 = process_time()
		res1 = left_to_right_square_and_multiply(a, b, p, n)
		sum1 += process_time() - c1
		c1 = process_time()
		res2 = left_to_right_square_and_multiply_always(a, b, p, n)
		sum2 += process_time() - c1
		c1 = process_time()
		res3 = montgomery_ladder(a, b, p, n)
		sum3 += process_time() - c1
		c1 = process_time()
		res4 = pow(a, b, p)
		sum4 += process_time() - c1
		if res1 != res2 != res3 != res4:
			print("res1:", res1)
			print("res2:", res2)
			print("res3:", res3)
			print("res4:", res4)
			print()
	print("res:\t", sum1, "\t", sum2, "\t", sum3, "\t", sum4)

	print("\nModular Inverse")
	sum1 = 0
	sum2 = 0
	sum3 = 0
	for _ in range(100):
		n = randrange(64, 128)
		p = randrange(1<<(n-1), 1<<n)
		p = next_prime(p)
		a = randrange(2, 1<<n)
		b = randrange(2, 1<<n)
		u, v, d = extended_euclidan(a, b)
		if a * u + b * v != d:
			print(False)
		else:
			c1 = process_time()
			res1 = pow(a, p - 2, p)
			sum1 += process_time() - c1
			c1 = process_time()
			res2, soak, soak = extended_euclidan(a, p)
			sum2 += process_time() - c1
			c1 = process_time()
			res3 = pow(a, -1, p)
			sum3 += process_time() - c1
			if res1 != (res2 % p) != res3:
				print("res1:", res1)
				print("res2:", res2)
				print("res3:", res3)
				print()
	print("res:\t", sum1, "\t", sum2, "\t", sum3)

	# From here on MNS
	print("\nMNS")
	l, n, k, lam, ksi = 272, 5, 60, 2, [0, 0, 0, 0, -1]
	p, gamma, rho, M = create_mns(l, n, k, lam, ksi)
	mns = (p, n, gamma, rho)
	print(mns)
	for _ in range(1):
		a = randrange(p)
		A = naive_convert_to_mns(a, *mns)
		if horner_modulo(A, gamma, p) != a:
			print(False)
			print(horner_modulo(A, gamma, p) - a)
		b = randrange(p)
		B = naive_convert_to_mns(b, *mns)
		A = [3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128]
		B = [1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068]
		a, b = horner_modulo(A, gamma, p), horner_modulo(B, gamma, p)
		Bcopy = B
		C = mns_mod_mult(A, B, *mns, lam)
		c = (a * b) % p
		if horner_modulo(C, gamma, p) != c:
			print("Error:")
			print("c")
			print(c)
			print("c_check")
			print(horner_modulo(C, gamma, p))
			print("difference")
			print(c - horner_modulo(C, gamma, p))
			break

	# Babai rounding method for coefficient reduction
	B = [[-pow(gamma, i, p) if j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]
	B[0][0] = p
	Betoile = [[int(sum([B[i][j] * B[k][j] for j in range(n)]) % p) for k in range(n)] for i in range(n)]
	Betoile[0][0] = p
	Cprime = babai_coefficient_reduction(C, *mns, B, Betoile)
	if horner_modulo(C, gamma, p) != horner_modulo(Cprime, gamma, p):
		print("WRONG")
		print(horner_modulo(C, gamma, p))
		print(horner_modulo(Cprime, gamma, p))
		print(int(horner_modulo(Cprime, gamma, p)) - int(horner_modulo(C, gamma, p)))
		print(gamma)
	print("Babai rounding causes smaller coefficients:")
	CTEST = [abs(Cprime[i]) <= abs(C[i]) for i in range(n)]
	print(CTEST)
	B = IntegerLattice(matrix(ZZ, B)).LLL()
	Betoile = [[int(sum([B[i][j] * B[k][j] for j in range(n)]) % p) for k in range(n)] for i in range(n)]
	Cseconde = babai_coefficient_reduction(C, *mns, B, Betoile)
	print("Babai rounding causes smaller coefficients after LLL:")
	CTEST = [abs(Cseconde[i]) <= abs(C[i]) for i in range(n)]
	print(CTEST)
	if horner_modulo(C, gamma, p) != horner_modulo(Cseconde, gamma, p):
		print("WRONG AGAIN")

	# Barrett like
	# TODO: if I have time later on (probably not)

	# Montgomery
	print("Montgomery reduction")
	B = list(B)
	RingPoly = PolynomialRing(ZZ, 'X')
	E = RingPoly("X^5 - 2")
	phi = 2**64
	for lig in B:
		for elem in lig:
			if abs(elem) > rho:
				print("Something wrong with LLL'd matrix")
		M = RingPoly(list_to_poly(lig))
		val, M1, soak = xgcd(M, E)
		val = ZZ(val)
		if val & 1:  # if val is even it won't have a gcd of 1 with phi.
			M1 = (M1 * ZZ(pow(val, -1, phi)) % phi)
			print(((M * M1) % E) % phi)
			print(M(gamma) % p == 0)
			break
	M = list(M)
	print("M")
	print(M)
	M1 = list(M1)
	print("M1")
	M1 = [(int(M1[i]) * -1) % phi for i in range(n)]
	M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]
	print(M1)
	print("Ctierce")
	Ctierce = montgomery_like_coefficient_reduction(C, *mns, lam, phi, M, M1)
	print(Ctierce)
	print(horner_modulo(Ctierce, gamma, p) * phi % p == horner_modulo(C, gamma, p))
	print("Montgomery-like causes smaller coefficients:")
	CTEST = [abs(Ctierce[i]) <= abs(C[i]) for i in range(n)]
	print(CTEST)
	print("Montgomery-like causes coefficients under rho:")
	CTEST = [abs(Ctierce[i]) <= rho for i in range(n)]
	print(CTEST)
	print("A")
	print(A)
	print("B")
	print(Bcopy)
	print("C")
	print(Ctierce)
	print(horner_modulo([-45427161615311394, 15593644179898433, 7040936443281178, 18209456034359938, -2124383526193652], gamma, p) == horner_modulo(Ctierce, gamma, p))
	print(amns_montg_mult(A, Bcopy, p, n, gamma, rho, lam, phi, M, M1))
	print("MM1")
	print([elem * lam for elem in M])
	print([elem * lam for elem in M1])

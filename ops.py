from random import randrange
from time import process_time
from sage.all import next_prime, factor, matrix
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

def create_mns(n, k, l, lam=2):
	while True:
		ksi = [0] * n
#		for _ in range(randrange(1, 3)):
#			ksi[randrange(n)] = randrange(-1, 2, 2)
		ksi[4] = -1
		M = [[ksi[(i - j) % (n)] * ((i - j < 0) * lam + (i - j >= 0)) for i in range(n)] for j in range(n)]
		M = [[-M[i][j] if i != j else (1<<k) -M[i][i] for i in range(n)] for j in range(n)]
		det = gauss_pivot(M)
#		det = matrix(M).determinant()
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

def mns_mod_mult(A, B, n, lam):
	R = [0] * n
	for i in range(n):
		for j in range(1, n - i):
			print("lambda", i+j, n-j)
			R[i] += A[i + j] * B[n - j]
		R[i] *= lam
		for j in range(i + 1):
			print(j, i-j)
			R[i] += A[j] * B[i-j]
	return R


def external_reduction(C, n, lam):
	R = [0] * n
	for i in range(n - 1):
		R[i] = C[i] + lam * C[n + i]
	R[n - 1] = C[n - 1]
	return R

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
	for _ in range(100):
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

	print("\nMNS")
	print(create_mns(5, 60, 272))

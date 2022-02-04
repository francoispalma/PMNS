from random import randrange
from time import process_time

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

def pseudo_mersenne_modular_mult(a, b, p, n):
	c = 

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

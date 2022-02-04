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

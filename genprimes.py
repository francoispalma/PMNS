from sage.all import previous_prime

with open("primes4096", "w+") as f:
	f.write("PRIMES4096 = []")
	p = 2**4096
	for _ in range(1000):
		p = previous_prime(p)
		f.write("\nPRIMES4096 += [" + str(p) + "]")

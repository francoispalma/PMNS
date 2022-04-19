with open("primes2048", "w+") as f:
	f.write("PRIMES2048 = []")
	p = 2**2048
	for _ in range(1000):
		p = previous_prime(p
		f.write("\nPRIMES2048 += [" + str(p) + "]")

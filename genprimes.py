from sage.all import previous_prime

with open("generated/primes8192.py", "w+") as f:
	f.write("PRIMES8192 = []")
	p = 2**8192
	for _ in range(1000):
		p = previous_prime(p)
		f.write("\nPRIMES8192 += [" + str(p) + "]")

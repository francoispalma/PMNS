import sys

from sage.all import random_prime


if __name__ == "__main__":
	if len(sys.argv) >= 2:
		strg = sys.argv[1]
		try:
			psize = int(strg)
			if psize < 32:
				print("Prime Size not handled")
				exit()
			with open(f"generated/primes{strg}.py", "w+") as f:
				f.write(f"PRIMES{strg} = []")
				primeset = set()
				while len(primeset) < 1000:
					p = random_prime(2**psize - 1, lbound= 2**(psize - 1))
					if p not in primeset:
						primeset.add(p)
						f.write(f"\nPRIMES{strg} += [" + str(p) + "]")
		except ValueError:
			print("Invalid arguments")
	else:
		print("Not enough arguments: Psize")

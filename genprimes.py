import sys

from sage.all import previous_prime


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
				p = 2**psize
				for _ in range(1000):
					p = previous_prime(p)
					f.write(f"\nPRIMES{strg} += [" + str(p) + "]")
		except ValueError:
			print("Invalid arguments")
	else:
		print("Not enough arguments: Psize")

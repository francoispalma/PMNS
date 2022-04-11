from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor
from math import ceil

from precalcs128 import do_precalcs as do_precalcs128
from precalcs import do_precalcs
from primes1024 import PRIMES1024

if __name__ == "__main__":
	print("pmnsdict = {}")
	for i in range(1000):
		p = PRIMES1024[i]
		K = GF(p)
		polK = PolynomialRing(K, 'X')
		n = 17
		flag = False
		while True:
			for lam in range(2, 8):
				E = polK("X^" + str(n) +" - " + str(lam))
				Eprime = polK("X^" + str(n) +" + " + str(lam))
				fs = factor(E)
				fsprime = factor(Eprime)
				if fs[0][0].degree() == 1:
					flag = True
					lamb = lam
					break
				if fsprime[0][0].degree() == 1:
					flag = True
					fs = fsprime
					lamb = -lam
					break
			if flag == True:
				break
			n += 2
		print("pmnsdict[" + str(p) + "] = (" + str(p) + ", " + str(n) + ", " + str(fs[0][0][0]) + ", " + str(lamb) + ")")

from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor
from math import ceil

from precalcs128 import do_precalcs as do_precalcs128
from precalcs import do_precalcs
from primes1024 import PRIMES1024

#generate P
#K = GF(P)
#polK.<X> = K[]
#P = X^n - lambda
#factor(P)
#lone factor => nth root of lambda.

if __name__ == "__main__":
#	p = 2**1024
#	print("[", end="")
#	for i in range(999):
#		p = previous_prime(p)
#		print(p, end=",\n")
#	print(previous_prime(p), end="")
#	print("]")
	print("pmnsdict = {}")
	print("pmns128dict = {}")
	for i in range(1000):
		p = PRIMES1024[i]
		K = GF(p)
		polK = PolynomialRing(K, 'X')
		n = 9
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
			n += 1
		#assert pow(fs[0][0][0], n, p) == (p - lamb) % p
		print("pmns128dict[" + str(p) + "] = (" + str(p) + ", " + str(n) + ", " + str(fs[0][0][0]) + ", " + str(lamb) + ")")
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
			n += 1
		#assert pow(fs[0][0][0], n, p) == (p - lamb) % p
		print("pmnsdict[" + str(p) + "] = (" + str(p) + ", " + str(n) + ", " + str(fs[0][0][0]) + ", " + str(lamb) + ")")

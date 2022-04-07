from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor
from math import ceil

from precalcs128 import do_precalcs as do_precalcs128
from precalcs import do_precalcs

#generate P
#K = GF(P)
#polK.<X> = K[]
#P = X^n - lambda
#factor(P)
#lone factor => nth root of lambda.

if __name__ == "__main__":
	p = 2**1024
	print("[", end="")
	for i in range(999):
		p = previous_prime(p)
		print(p, end=",\n")
	print(previous_prime(p), end="")
	print("]")
#		K = GF(p)
#		polK = PolynomialRing(K, 'X')
#		E = polK("X^9 - 2")
#		fs = factor(E)
#		print(fs[0][0].degree())


from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

from primes2048 import PRIMES2048

if __name__ == "__main__":
	print("pmnsdict = {}")
	for i in range(1000):
		p = PRIMES2048[i]
		K = GF(p)
		polK = PolynomialRing(K, 'X')
		n = 33
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
				flag = False
				gamma = fs[0][0][0]
				B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
				B = list(IntegerLattice(matrix(ZZ, B)).LLL())
				w = 1 + (n - 1) * abs(lam)
				__tmp = int(2 * w * max(max(B)))
				rho = ceil(__tmp.bit_length())
				if rho <= 64:
					break
			n += 2
		print("pmnsdict[" + str(p) + "] = (" + str(p) + ", " + str(n) + ", " + str(fs[0][0][0]) + ", " + str(lamb) + ")")

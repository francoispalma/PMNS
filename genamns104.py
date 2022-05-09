from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

from primes4096 import PRIMES4096 as PRIMES

# ||M||inf <= phi/((lambda*n*2)**2)
# min(||M||inf) = P^1/n

if __name__ == "__main__":
	print("pmns104dict = {}")
	power = 4096
	phi = 104
	PHI = 2**phi
	init_n = (power // phi) | 1
	while 2**(power/init_n) >= PHI/((2*init_n*2)**2):
		init_n += 2
	for i in range(1000):
		p = PRIMES[i]
		K = GF(p)
		polK = PolynomialRing(K, 'X')
		n = init_n
		flag = False
		while True:
			POWERN = 2**(power/n)
			for lam in range(2, 8):
				w = 1 + (n - 1) * lam
				if PHI <= POWERN * 2 * w:
					break
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
				__tmp = int(2 * w * max(max(B)))
				rho = ceil(__tmp.bit_length())
				if rho < 2 * w * phi:
					break
			n += 2
		print("pmns104dict[" + str(p) + "] = (" + str(p) + ", " + str(n) + ", " + str(fs[0][0][0]) + ", " + str(lamb) + ")")

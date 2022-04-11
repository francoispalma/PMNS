from sage.all import matrix, ZZ, PolynomialRing
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

from generatedpmns import pmnsdict, pmns128dict
from primes1024 import PRIMES1024 as primes
from precalcs128 import do_precalcs as do_precalcs128
from precalcs import do_precalcs

for i in range(1000):
	(p, n, gamma, lam) = pmns128dict[primes[i]]
	#do_precalcs(p, n, gamma, lam)
	B = [[p if (i, j) == (0, 0) else -pow(gamma, i, p) if i != 0 and j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]
	B = list(IntegerLattice(matrix(ZZ, B)).LLL())
	RingPoly = PolynomialRing(ZZ, 'X')
	E = RingPoly("X^" + str(n) + " - (" + str(lam) + ")")
	w = 1 + (n - 1) * abs(lam)
	__tmp = int(2 * w * max(max(B)))
	rho = ceil(__tmp.bit_length())
	if rho > 128:
		print(i, rho)

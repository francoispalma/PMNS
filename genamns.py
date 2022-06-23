from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ
from sage.modules.free_module_integer import IntegerLattice
from math import ceil, log2
import sys

from commonpmns import primesdict, handledphis

# ||M||inf <= phi/((lambda*n*2)**2)
# min(||M||inf) = P^1/n

def gen_amns(power, sphi, start=0):
	primes = primesdict[power]
	if start == 0:
		print("pmns" + sphi + "dict = {}")
	if sphi == "":
		phi = 64
	else:
		phi = int(sphi)
	PHI = 2**phi
	init_n = (power // phi) | 1
	while 2**(power/init_n) >= PHI/((2*init_n*2)**2):
		init_n += 2
	for i in range(len(primes)):
		p = primes[i]
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
				__tmp = int(2 * w * max([max(Line) for Line in B])))
				rho = ceil(__tmp.bit_length())
				if rho < phi - 1 - log2(w):
					break
			n += 2
		print(f"pmns{sphi}dict[{p}] = ({p}, {n}, {fs[0][0][0]}, {lamb}, {rho}, {B})")

if __name__ == "__main__":
	if len(sys.argv) >= 2:
		try:
			psize = int(sys.argv[1])
			if psize not in list(primesdict.keys()):
				print("Prime Size not handled")
				exit()
			try:
				phi = sys.argv[2]
			except IndexError:
				phi = ""
			if phi not in handledphis:
				print("Value of Phi not handled")
				exit()
			elif phi == "64":
				phi = ""
			if len(sys.argv) >= 4:
				start = int(sys.argv[3])
			else:
				start = 0
			gen_amns(psize, phi, start)
		except ValueError:
			print("Invalid arguments")
	else:
		print("Not enough arguments: Psize [Phi] [start]")

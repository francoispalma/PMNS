from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import ceil, log2
import sys

from ops import list_to_poly
from commonpmns import primesdict, handledphis

# ||M||inf <= phi/((lambda*n*2)**2)
# min(||M||inf) = P^1/n

def gen_amns(power, sphi, start=0):
	primes = primesdict[power]
	if start == 0:
		print("pmns" + sphi + "dict = {}")
	else:
		print()
	if sphi == "":
		phi = 64
	else:
		phi = int(sphi)
	PHI = 2**phi
	init_n = (power // phi) | 1
	while 2**(power/init_n) >= PHI/((2*init_n*2)**2):
		init_n += 2
	for i in range(start, len(primes)):
		p = primes[i]
		K = GF(p)
		polK = PolynomialRing(K, 'X')
		n = init_n
		flag = False
		while True:
			POWERN = 2**(power/n)
			# We try for each value of lambda
			for lam in range(2, 8):
				w = 1 + (n - 1) * lam
				if PHI <= POWERN * 2 * w:
					break

				# We define our external reduction polynomial
				E = polK("X^" + str(n) +" + " + str(lam))
				fs = factor(E)
				if fs[0][0].degree() == 1:  # if the degree is one, we have a solution
					flag = True
					lamb = lam
					break

				# We also have to try with negative values
				Eprime = polK("X^" + str(n) +" - " + str(lam))
				fsprime = factor(Eprime)
				if fsprime[0][0].degree() == 1:
					flag = True
					fs = fsprime
					E = Eprime
					lamb = -lam
					break
			if flag == True:
				flag = False
				gamma = fs[0][0][0]

				RingPoly = PolynomialRing(ZZ, 'X')
				E = RingPoly("X^" + str(n) + " - (" + str(lam) + ")")

				# We calculate the base matrix
				B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
				B = list(IntegerLattice(matrix(ZZ, B)).LLL())

				# Then we calculate rho
				__tmp = int(2 * w * max([max(Line) for Line in B]))
				rho = ceil(__tmp.bit_length())
				if rho < phi - 1 - log2(w):
					# We then try to find a valid M
					for lig in B:
						# We set M as one line of the base matrix
						M = RingPoly(list_to_poly(lig))
						# We then get the d and u of the au + bv = d from extended euclid
						val, M1, soak = xgcd(M, E)
						# We make sure to get out of any sage specific field
						val = ZZ(val)

						if val & 1:  # if val is even it won't have a gcd of 1 with phi
							M1 = (M1 * ZZ(pow(val, -1, PHI)) % PHI)
							# We check that M1 is properly constructed, if it is we don't look further
							if ((M * M1) % E) % PHI == 1 and M(gamma) % p == 0:
								break

					# We convert out of the polynomial Ring
					M = list(M)
					M1 = list(M1)

					# We switch M1 from M^-1 to -M^-1
					M1 = [(int(M1[i]) * -1) % PHI for i in range(n)]

					# We then reduce M1's coefficients
					M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else M1[i] for i in range(n)]
					break
			n += 2
		print(f"pmns{sphi}dict[{p}] = ({p}, {n}, {gamma}, {lamb}, {rho}, {M}, {M1})")

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

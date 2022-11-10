from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import log2
import sys

from commonpmns import primesdict, handledphis


def list_to_poly(L: list) -> str:
	"""Utility function to transform lists to strings to convert to polynomials
	in the format sage likes.
	L: The list that represents our polynomials.
	"""
	rstr = str(L[0])
	for i in range(1, len(L)):
		rstr += " + (" + str(L[i]) + ") * X" + ("^" + str(i) if i != 1 else "")
	return rstr

def infinite_norm_of_matrix(Matrix):
	return max([sum([abs(elem) for elem in Line]) for Line in Matrix])

def norm_one_of_matrix(Matrix):
	return max([sum([abs(Matrix[i][j]) for i in range(len(Matrix))]) for j in range(len(Matrix))])

def maingen(psize, sphi, start=0, polyv=True, delta=0):
	primes = primesdict[psize]
	if start == 0:
		print("pmns" + sphi + "dict = {}")
	else:
		print()
	if sphi == "":
		phi = 64
	else:
		phi = int(sphi)
	for i in range(start, len(primes)):
		p = primes[i]
		p, n, gamma, lamb, rho, M_or_B, M1_or_B1 = gen_amns(p, phi, polyv, delta)
		print(f"pmns{sphi}dict[{p}] = ({p}, {n}, {gamma}, {lamb}, {rho}, {M_or_B}, {M1_or_B1})")

def gen_amns(p, phi=64, polyv=True, delta=0):
	PHI = 2**phi
	psize = int(log2(p))
	# We try to find the minimal n for a given p

	# This is the smallest theoretical value
	n = (psize // phi)

	#PHI >= 4wn1B >= 4wp^(1/n)
	#PHI >= 4w2**(log2(p)/n)
	#PHI >= 4(2n - 1)2**(log2(p)/n)
	#phi >= 2 + log2(2n-1) + log2(p)/n
#	while phi < 2 + log2(2*n-1) + log2(p)/n:
	while PHI < 4 * (2*n - 1) * (2**(log2(p)/n)) * ((delta + 1)**2):
		n += 1
	K = GF(p)
	polK = PolynomialRing(K, 'X')
	flag = False
	while True:
		POWERN = 2**(psize/n)
		# We try for each value of lambda
		for lam in range(2, 8):
			w = 1 + (n - 1) * lam
			if PHI <= POWERN * 2 * w:
				break

			if (n & 1):
				d, u, v = xgcd(n, p-1)
				if d == 1:
					flag = True
					potgamma = [pow(2, u, p)]
					break

			# We define our external reduction polynomial X^n - lambda
			E = polK("X^" + str(n) +" - " + str(lam))
			fs = factor(E)
			i = 0
			potgamma = []
			while i < len(fs) and fs[i][0].degree() == 1:
				potgamma += [-fs[i][0][0] % p]
				flag = flag or (abs(int(((pow(potgamma[-1], n, p) + 10) % p) - 10)) == lam)
				i += 1
			if flag:
				break

			# We also have to try with negative values so X^n + lambda in effect
			E = polK("X^" + str(n) +" + " + str(lam))
			fs = factor(E)
			i = 0
			potgamma = []
			while i < len(fs) and fs[i][0].degree() == 1:
				potgamma += [-fs[i][0][0] % p]
				flag = flag or (abs(int(((pow(potgamma[-1], n, p) + 10) % p) - 10)) == lam)
				i += 1
			if flag:
				break

		if flag == True:
			flag = False
			# We do this so that we get negative values instead of p sized lambdas
			lamb = int(((pow(potgamma[0], n, p) + 10) % p)) - 10

			RingPoly = PolynomialRing(ZZ, 'X')
			E = RingPoly("X^" + str(n) + " - (" + str(lamb) + ")")

			for gamma in potgamma:
				# We calculate the base matrix
				B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
				B = list((matrix(ZZ, B).LLL()))

				# We want PHI => 4wnorm1(B) because we want PHI => 2wrho
				# and rho => 2norm1(B)
				n1B = norm_one_of_matrix(B)
				if PHI >= 4 * w * n1B * ((delta + 1)**2):

					if not polyv:
						B1 = list(matrix(B).inverse() % PHI)
						__tmp = int(2 * n1B)
						rho = __tmp.bit_length()  # We take a bigger rho for now
						return p, n, gamma, lamb, rho, B, B1

					#TODO: delete it, this is the old method.
					elif polyv:
						# We then try to find a valid M
						for lig in B:
							# We set M as one line of the base matrix
							M = RingPoly(list_to_poly(lig))
							# We then get the d and u of the au + bv = d from extended euclid
							val, M1, soak = xgcd(M, E)
							# We make sure to get out of any sage specific field
							try:
								val = ZZ(val)
							except TypeError:
								val = 0

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
						__tmp = int(2 * w * max([abs(elem) for elem in M]))
						rho = __tmp.bit_length()
						if rho <= phi - 1 - log2(w):
							return p, n, gamma, lamb, rho, M, M1

				
		n += 1
	

if __name__ == "__main__":
	args = sys.argv.copy()
	opts = [arg for arg in args if arg.startswith("-")]
	if "-d" not in opts and "--delta" not in opts:
		delta = 0
	else:
		if "-d" in opts:
			op = "-d"
		else:
			op = "--delta"
		try:
			delta = int(args[args.index(op) + 1])
		except IndexError:
			print("Error in use of option:", op)
			print("Not enough arguments")
			exit()
		except ValueError:
			print("Invalid value for delta:", args[args.index(op) + 1])
			exit()
		args.pop(args.index(op) + 1)
	assert delta >= 0, f"Invalid value for delta: {delta}"
	arguments = [arg for arg in args if not arg.startswith("-")]
	if len(arguments) >= 2:
		try:
			psize = int(arguments[1])
		except ValueError:
			print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
			exit()
		if psize not in list(primesdict.keys()):
			print("Prime Size not handled")
			exit()
		try:
			phi = arguments[2]
		except IndexError:
			phi = ""
		if phi not in handledphis:
			print("Value of Phi not handled:", phi)
			exit()
		elif phi == "64":
			phi = ""
		if len(arguments) >= 4:
			start = int(arguments[3])
		else:
			start = 0
		if "-B" in opts or "--Base" in opts:
			maingen(psize, phi, start, False, delta)
		else:
			maingen(psize, phi, start, True, delta)
	else:
		print("Not enough arguments: Psize [Phi] [start]\n\nOPTIONS")
		print("\t--Base, -B: generate base reduction version instead of polynomial.")
		print("\t--delta, -d: use specified delta for number of free additions.")

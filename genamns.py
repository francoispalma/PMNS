from sage.all import next_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import ceil, log2
import sys

from ops import list_to_poly
from commonpmns import primesdict, handledphis

# ||M||inf <= phi/((lambda*n*2)**2)
# min(||M||inf) = P^1/n

def infinite_norm_of_matrix(Matrix):
	return max([sum([abs(elem) for elem in Line]) for Line in Matrix])

def norm_one_of_matrix(Matrix):
	return max([sum([abs(Matrix[i][j]) for i in range(len(Matrix))]) for j in range(len(Matrix))])

def maingen(power, sphi, start=0, polyv=True):
	primes = primesdict[power]
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
		p, n, gamma, lamb, rho, M_or_B, M1_or_B1 = gen_amns(p, phi, polyv)
		print(f"pmns{sphi}dict[{p}] = ({p}, {n}, {gamma}, {lamb}, {rho}, {M_or_B}, {M1_or_B1})")

def gen_amns(p, phi=64, polyv=True):
	PHI = 2**phi
	power = int(log2(p))
	n = (power // phi)
	while power >= n * phi - n * (4 + 2*log2(n)):
		n += 1
	K = GF(p)
	polK = PolynomialRing(K, 'X')
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
			try:
				fs = factor(E)
				if fs[0][0].degree() == 1:  # if the degree is one, we have a solution
					lamb = lam
					if(pow(fs[0][0][0], n, p) == lamb):
						flag = True
						break
			except:
				pass

			# We also have to try with negative values
			Eprime = polK("X^" + str(n) +" - " + str(lam))
			fsprime = factor(Eprime)
			if fsprime[0][0].degree() == 1:
				fs = fsprime
				E = Eprime
				lamb = -lam
				if(pow(fs[0][0][0], n, p) == lamb):
					flag = True
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
			#__tmp = int(2 * w * infinite_norm_of_matrix(B))
			if polyv:
#				__tmp = int(2 * w * max([abs(B[i][j]) for i in range(n) for j in range(n)]))
				__tmp = int(2 * w * min([max(lig) for lig in B]))
				rho = ceil(__tmp.bit_length())
				if rho < phi - 1 - log2(w):
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
					rho = ceil(__tmp.bit_length())
					return p, n, gamma, lamb, rho, M, M1

			else:
				__tmp = int(2 * norm_one_of_matrix(B))
				rho = ceil(__tmp.bit_length())
				if rho < phi - 1 - log2(w):
					B1 = list(matrix(B).inverse() % PHI)
					return p, n, gamma, lamb, rho, B, B1

		n += 1
	

if __name__ == "__main__":
	arguments = sys.argv.copy()
	opts = [arg for arg in arguments if arg.startswith("-")]
	arguments = [arg for arg in arguments if not arg.startswith("-")]
	if len(arguments) >= 2:
		try:
			psize = int(arguments[1])
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
				maingen(psize, phi, start, False)
			else:
				maingen(psize, phi, start, True)
		except ValueError:
			print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
	else:
		print("Not enough arguments: Psize [Phi] [start]\n\nOPTIONS\n\t--Base, -B: generate base reduction version instead of polynomial.")

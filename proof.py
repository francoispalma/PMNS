from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import ceil
from random import randrange
from ast import literal_eval

from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from ops import list_to_poly, horner_modulo, amns_montg_mult
from pyparams import p, n, gamma, lam, phi

#amns = (p, n, gamma, rho, lam) = (6152896135288560374679945371974689688835168151742564408104565373600581564260451457, 5, 220855883097298041197912187592864814478435487109452369765200775161577472, 2305843009213693952, 2)
#phi = 2**64
#M = [9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940]
#M1 = [7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353]

if __name__ == "__main__":
	# We calculate the base matrix
	B = [[p if (i, j) == (0, 0) else -pow(gamma, i, p) if i != 0 and j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]

	# We apply LLL to it
	B = list(IntegerLattice(matrix(ZZ, B)).LLL())

	# We define our external reduction polynomial
	RingPoly = PolynomialRing(ZZ, 'X')
	E = RingPoly("X^" + str(n) + " - (" + str(lam) + ")")

	# We find our w factor to get rho
	w = 1 + (n - 1) * abs(lam)
	__tmp = int(2 * w * max([max(Line) for Line in B])))
	RHO = ceil(__tmp.bit_length())
	rho = 2**RHO
	
	# We then try to find a valid M
	for lig in B:
		# We check if something went wrong
		for elem in lig:
			if abs(elem) > rho:
				print("Something wrong with LLL'd matrix")
				exit()

		# We set M as one line of the base matrix
		M = RingPoly(list_to_poly(lig))
		# We then get the d and u of the au + bv = d from extended euclid
		val, M1, soak = xgcd(M, E)
		# We make sure to get out of any sage specific field
		val = ZZ(val)

		if val & 1:  # if val is even it won't have a gcd of 1 with phi
			M1 = (M1 * ZZ(pow(val, -1, phi)) % phi)
			# We check that M1 is properly constructed, if it is we don't look further
			if ((M * M1) % E) % phi == 1 and M(gamma) % p == 0:
				break

	# We convert out of the polynomial Ring
	M = list(M)
	M = M + (n - len(M)) * [0]
	M1 = list(M1)
	M1 = M1 + (n - len(M1)) * [0]

	# We switch M1 from M^-1 to -M^-1
	M1 = [(int(M1[i]) * -1) % phi for i in range(n)]

	# We then reduce M1's coefficients
	M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]
	
	counter = 0
	print("Starting")
	with open("log", "r") as f:
		while True:
			try:
				a = literal_eval(f.readline()[:-1])
				a = [elem - phi if elem > (phi >> 1) else elem for elem in a]
				b = literal_eval(f.readline()[:-1])
				b = [elem - phi if elem > (phi >> 1) else elem for elem in b]
				c = literal_eval(f.readline()[:-1])
				c = [elem - phi if elem > (phi >> 1) else elem for elem in c]
				c_check = amns_montg_mult(a, b, p, n, gamma, rho, lam, phi, M, M1)
				if horner_modulo(c, gamma, p) != horner_modulo(c_check, gamma, p):
					counter += 1
					#print("False")
				for elem in c:
					if abs(elem) >= rho:
						print("More than Rho")
			except SyntaxError:
				break
	print("Finished")
	print("counter:", counter)

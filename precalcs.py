from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

from ops import list_to_poly, montgomery_like_coefficient_reduction, horner_modulo, amns_montg_mult
from convert import montgomery_convert_to_mns, rho_div_convert_to_mns
from proof import p, n, gamma, lam, phi

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, phi):
	print("#ifndef PMNS_PARAMS_H_INCLUDED\n#define PMNS_PARAMS_H_INCLUDED\n")

	# We calculate the base matrix
	B = [[p if (i, j) == (0, 0) else -pow(gamma, i, p) if i != 0 and j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]

	# We apply LLL to it
	B = list(IntegerLattice(matrix(ZZ, B)).LLL())

	# We define our external reduction polynomial
	RingPoly = PolynomialRing(ZZ, 'X')
	E = RingPoly("X^" + str(n) + " - (" + str(lam) + ")")

	# We find our w factor to get rho
	w = 1 + (n - 1) * abs(lam)
	__tmp = int(2 * w * max(max(B)))
	rho = ceil(__tmp.bit_length())
	print("#define RHO", rho)
	rho = 2**rho

	print("#define N", str(n) + "\n#define LAMBDA", str(lam) +"\n")
	print("// Various values you can get from precalc")

	# We then try to find a valid M
	for lig in B:
		# We check if something went wrong
		for elem in lig:
			if abs(elem) > rho:
				print("Something wrong with LLL'd matrix")
				return

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
	M1 = list(M1)

	# We switch M1 from M^-1 to -M^-1
	M1 = [(int(M1[i]) * -1) % phi for i in range(n)]

	# We then reduce M1's coefficients
	M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]

	# We print
	print("const int64_t M[] = {" + str(M)[1:-1] + "},")
	print("\tM1[] = {" + str(M1)[1:-1] + "},")
	print("\tMLambda[] = {" + str([elem * lam for elem in M])[1:-1] + "},")
	M1L = [elem * lam for elem in M1]
	M1L = [M1L[i] - phi if M1L[i] >= (phi >> 1) else M1L[i] for i in range(n)]
	print("\tM1Lambda[] = {" + str(M1L)[1:-1] + "},")

	# We then get the Pi and print it
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), p, n, gamma, rho, lam, phi, M, M1, phinmoinsun) for i in range(n)]
	print("\t__Pi__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str(Pi[i])[1:-1] + "},")
	print("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n")

	# Powers of gamma next
	string = "G"
	print("const _poly ", end="")
	g = gamma
	for i in range(1, n):
		string = string + str(i) + ", G"
		tmp = convert_to_int_tabs(int(g))
		if i != 1:
			print("\t", end="")
		print("G" + str(i) + " = { .deg = " + str(len(tmp)) + ",")
		print("\t\t.t = (int64_t[]) {", end="")
		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
		print(tmp + "} }", end="")
		if i != n - 1:
			print(",")
		else:
			print(";")
		g = g * gamma % p
	
	print("_poly __P__ = { .deg = " + str(n) + ",")
	tmp = convert_to_int_tabs(p)
	tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
	print("\t\t.t = (int64_t[]) {" + tmp + "} },")
	print("\tGi[] = {" + string[:-3] + "};\n")

	print("#endif")

if __name__ == "__main__":
#	p = 5524969863260095610495186344939419738027465949850580467176395575832917506871951737255621939342449907372936940924410124929668828406973361712220148691590192943
#	n = 12
#	gamma = 4636652768285458569062878611317602841499088282530798143147640460012527948722544350531006403703239417663454147165839392083945105225074778235517366830265723868
#	lam = 3
#	phi = 2**64
	do_precalcs(p, n, gamma, lam, phi)

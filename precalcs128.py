from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

from ops import list_to_poly, montgomery_like_coefficient_reduction, horner_modulo, amns_montg_mult
from convert import montgomery_convert_to_mns, rho_div_convert_to_mns
#from proof import p, n, gamma, lam, phi
from pyparams128 import p, n, gamma, lam

phi = 2**128

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam):
	print("#ifndef PMNS_PARAMS128_H_INCLUDED\n#define PMNS_PARAMS128_H_INCLUDED\n")

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
	M = M + (n - len(M)) * [0]
	M1 = list(M1)
	M1 = M1 + (n - len(M1)) * [0]

	# We switch M1 from M^-1 to -M^-1
	M1 = [(int(M1[i]) * -1) % phi for i in range(n)]

	# We then reduce M1's coefficients
	M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]

	# We get M1Lambda
	M1L = [(elem * lam) % phi for elem in M1]
	M1L = [M1L[i] - phi if M1L[i] >= (phi >> 1) else M1L[i] for i in range(n)]

	# We cut up M, MLambda, M1 and M1lambda
	Mhi = [elem >> 64 for elem in M]
	Mlo = [elem % (2**64) for elem in M]
	MLambdalo = [(elem * lam) % 2**64 for elem in M]
	MLambdahi = [(elem * lam) >> 64 for elem in M]
	M1hi = [elem >> 64 for elem in M1]
	M1lo = [elem % (2**64) for elem in M1]
	M1Lambdalo = [elem % 2**64 for elem in M1L]
	M1Lambdahi = [elem >> 64 for elem in M1L]

	# We then get the Pi
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), p, n, gamma, rho, lam, phi, M, M1, phinmoinsun) for i in range(n)]

	# Then we print
	print("static const uint64_t Mlo[] = {" + str(Mlo)[1:-1].replace(",", "u,") + "u},")
	print("\tM1lo[] = {" + str(M1lo)[1:-1].replace(",", "u,") + "u},")
	print("\tMLambdalo[] = {" + str(MLambdalo)[1:-1].replace(",", "u,") + "u},")
	print("\tM1Lambdalo[] = {" + str(M1Lambdalo)[1:-1].replace(",", "u,") + "u},")
	print("\t__Pilo__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str([elem % 2**64 for elem in Pi[i]])[1:-1].replace(",", "u,") + "u},")
	print("\t\t{" + str([elem % 2**64 for elem in Pi[-1]])[1:-1].replace(",", "u,") + "u}\n\t};\n")

	print("static const int64_t Mhi[] = {" + str(Mhi)[1:-1] + "},")
	print("\tM1hi[] = {" + str(M1hi)[1:-1] + "},")
	print("\tMLambdahi[] = {" + str(MLambdahi)[1:-1] + "},")
	print("\tM1Lambdahi[] = {" + str(M1Lambdahi)[1:-1] + "},")
	print("\t__Pihi__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str([elem >> 64 for elem in Pi[i]])[1:-1] + "},")
	print("\t\t{" + str([elem >> 64 for elem in Pi[-1]])[1:-1] + "}\n\t};\n")

	# Powers of gamma next
	string = "G"
	print("static const _poly ", end="")
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

	print("static _poly __P__ = { .deg = " + str(n) + ",")
	tmp = convert_to_int_tabs(p)
	tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
	print("\t\t.t = (int64_t[]) {" + tmp + "} },")
	print("\tGi[] = {" + string[:-3] + "};\n")

	print("""
static inline void m1_mns128_mod_mult_ext_red_pre(unsigned __int128* restrict Rlo,
	const restrict poly128 A)
{
""")

	for i in range(n):
		print(f"Rlo[{i}] = (__int128)", end="")
		for j in range(1, n - i):
			print(f" ((__int128)A->lo[{i + j}] * {M1Lambdalo[n - j]}u) + " +
				f"((__int128) (LOW(A->lo[{i + j}] * {M1Lambdahi[n - j]}) + " +
				f"LOW(A->hi[{i + j}] * {M1Lambdalo[n - j]}u)) << 64)", end="")
			if j != n - i - 1:
				print(" +", end= "")
			else:
				print("+", end= "")

		for j in range(0, i + 1):
			print(f" ((__int128) A->lo[{j}] * {M1lo[i - j]}u) + ((__int128) " +
				f"(LOW(A->lo[{j}] * {M1hi[i - j]}) + LOW(A->hi[{j}] * " +
				f"{M1lo[i - j]}u)) << 64)", end="")
			if j != i:
				print(" +", end= "")

		print(";")

	print("}\n")

	print("""
static inline void m_mns128_mod_mult_ext_red_pre(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, const restrict poly128 A)
{
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux1, aux2, aux3;
""")

	for i in range(n):
		for j in range(1, n - i):
			print(f"\n\tA1B1 = ((__int128) A->hi[{i + j}] * {MLambdahi[n-j]});")
			print(f"\tA1B0 = ((__int128) A->hi[{i + j}] * {MLambdalo[n-j]}u);")
			print(f"\tA0B1 = ((__int128) A->lo[{i + j}] * {MLambdahi[n-j]});")
			print(f"\tA0B0 = ((__int128) A->lo[{i + j}] * {MLambdalo[n-j]}u);\n")
			print(f"""\taux3 = ((__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0));
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);
	aux1 = (__int128) HIGH(A1B1);
	
	tmplo = (__int128) LOW(A0B0) | (aux3 << 64);
	Rhi[{i}] += (__int128) aux2 + (aux1 << 64) +
		__builtin_add_overflow(Rlo[{i}], tmplo, Rlo + {i});""")

		for j in range(0, i + 1):
			print(f"\n\tA1B1 = ((__int128) A->hi[{j}] * {Mhi[i-j]});")
			print(f"\tA1B0 = ((__int128) A->hi[{j}] * {Mlo[i-j]}u);")
			print(f"\tA0B1 = ((__int128) A->lo[{j}] * {Mhi[i-j]});")
			print(f"\tA0B0 = ((__int128) A->lo[{j}] * {Mlo[i-j]}u);\n")
			print(f"""\taux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0) + LOW(A1B1);
	aux1 = (__int128) HIGH(A1B1);
	
	tmplo = (__int128) LOW(A0B0) | (aux3 << 64);
	Rhi[{i}] += (__int128) aux2 + (aux1 << 64) +
		__builtin_add_overflow(Rlo[{i}], tmplo, Rlo + {i});""")

	print("}\n")

	print("#endif")

if __name__ == "__main__":
#	p = 135235643069960614055763147653064061503447506836195743621526670176368184234720875133770858820476536434488752158042109408722131110950494765999584602783171012937442890496568352400298953673566640901275488915484314041362810104082501296367785554838617047203051373164870703338765011558032443485283660293658396373031
#	n = 9
#	gamma = 31655366034728588078624438514615301544129522139330836045091327852751402305472801538659044715201150659916554681038263189819329423658063171106562974285532277781484292697032498794895427121795477802587356170442668003146237461337340801098373462588444085927763403266750253334246732747899661111138239398839213291694
#	lam = 2
#	phi = 2**128
	do_precalcs(p, n, gamma, lam)

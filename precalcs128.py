import sys
from sage.all import matrix

from ops import list_to_poly, horner_modulo, amns_montg_mult as amns_montg_mult_poly, amns_montg_mult_base
from convert import montgomery_convert_to_mns as montgomery_convert_to_mns_poly, montgomery_convert_to_mns_base, rho_div_convert_to_mns as rho_div_convert_to_mns_poly, rho_div_convert_to_mns_base
from pyparams128 import p, n, gamma, lam, phi, rho, M_or_B, M1_or_B1
from commonpmns import pmnsdicts, primesdict

phi = 2**128

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, rho, M, M1):
	print("#ifndef PMNS_PARAMS128_H_INCLUDED\n#define PMNS_PARAMS128_H_INCLUDED\n")

	print("#define RHO", rho)
	rho = 2**rho

	print("#define N", str(n) + "\n#define LAMBDA", str(lam) +"\n")

	print("""
static inline _Bool add_overflow(unsigned __int128* restrict a, const unsigned __int128 b)
{
	//return __builtin_add_overflow(*a, b, a);
	const unsigned __int128 tmp = *a;
	*a += b;
	return *a < tmp;
}
""")
	
	print("// Various values you can get from precalc")

	auxstring = """\taux3 = (__int128) HIGH(A0B0) + LOW(A0B1) + LOW(A1B0);
	aux2 = (__int128) HIGH(aux3) + HIGH(A0B1) + HIGH(A1B0);
	
	tmplo = (__int128) LOW(A0B0) | (aux3 << 64);"""

	At = ["(__int128) (A->lo[", "(__int128) LOW(A->hi["]

	if type(M_or_B[0]) == int:
		print("#define M_or_B_is_M\n")
		M = M_or_B
		M1 = M1_or_B1

		montgomery_convert_to_mns = montgomery_convert_to_mns_poly
		amns_montg_mult = amns_montg_mult_poly
		rho_div_convert_to_mns = rho_div_convert_to_mns_poly

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

		# Then we print
		print("static const uint64_t Mlo[] = {" + str(Mlo)[1:-1].replace(",", "u,") + "u},")
		print("\tM1lo[N] = {" + str(M1lo)[1:-1].replace(",", "u,") + "u},")
		print("\tMLambdalo[N] = {" + str(MLambdalo)[1:-1].replace(",", "u,") + "u},")
		print("\tM1Lambdalo[N] = {" + str(M1Lambdalo)[1:-1].replace(",", "u,") + "u};")

		print("static const int64_t Mhi[N] = {" + str(Mhi)[1:-1] + "},")
		print("\tM1hi[N] = {" + str(M1hi)[1:-1] + "},")
		print("\tMLambdahi[N] = {" + str(MLambdahi)[1:-1] + "},")
		print("\tM1Lambdahi[N] = {" + str(M1Lambdahi)[1:-1] + "};")

		print("static const __int128 M1[N] = { ", end="")
		for i in range(n - 1):
			print(f"((__int128) LOW({M1hi[i]}) << 64) | {M1lo[i]}u, ", end="")
		print(f"((__int128) LOW({M1hi[n - 1]}) << 64) | {M1lo[n - 1]}u }},")
		print("\tM1Lambda[] = { ", end="")
		for i in range(n - 1):
			print(f"((__int128) LOW({M1Lambdahi[i]}) << 64) | {M1Lambdalo[i]}u, ", end="")
		print(f"((__int128) LOW({M1Lambdahi[n - 1]}) << 64) | {M1Lambdalo[n - 1]}u }};")

		# From here on, precalc functions, first with the multiplication by M1
		print("""
static inline void m1_or_b1_mns128_mod_mult_ext_red_pre(unsigned __int128* restrict Rlo,
	const restrict poly128 A)
{
""")

		for i in range(n):
			print(f"Rlo[{i}] = (__int128)", end="")
			for j in range(1, n - i):
				print(f" ((__int128) A->lo[{i + j}] * {M1Lambdalo[n - j]}u) + " +
					f"((__int128) (LOW(A->lo[{i + j}] * {M1Lambdahi[n - j]}) + " +
					f"LOW(A->hi[{i + j}] * {M1Lambdalo[n - j]}u)) << 64)", end="")
				print(" +", end= "")

			for j in range(0, i + 1):
				print(f" ((__int128) A->lo[{j}] * {M1lo[i - j]}u) + ((__int128) " +
					f"(LOW(A->lo[{j}] * {M1hi[i - j]}) + LOW(A->hi[{j}] * " +
					f"{M1lo[i - j]}u)) << 64)", end="")
				if j != i:
					print(" +", end= "")

			print(";")

		print("}\n")

		# Next is the multiplication by M
		print("""
static inline void m_or_b_mns128_mod_mult_ext_red_pre(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, const restrict poly128 A)
{
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux2, aux3;
""")

		Mt = [Mlo, Mhi]
		MLt = [MLambdalo, MLambdahi]

		for i in range(n):
			for j in range(1, n - i):
				for k in range(3, -1, -1):
					print(f"\tA{k//2}B{k&1} = {At[k//2]}{i + j}]) * {MLt[k&1][n-j]}" +
						f"{'' if k&1 else 'u'};")
				print(auxstring)
				print(f"""\tRhi[{i}] += (__int128) aux2 + A1B1 +
		add_overflow(Rlo + {i}, tmplo);\n""")

			for j in range(0, i + 1):
				for k in range(3, -1, -1):
					print(f"\tA{k//2}B{k&1} = {At[k//2]}{j}]) * {Mt[k&1][i-j]}" +
						f"{'' if k&1 else 'u'};")
				print(auxstring)
				print(f"""\tRhi[{i}] += (__int128) aux2 + A1B1 +
		add_overflow(Rlo + {i}, tmplo);\n""")

		print("}\n")

	elif type(M_or_B[0]) == tuple or type(M_or_B[0]) == list:
		print("#define M_or_B_is_B\n")
		B = M_or_B
		B1 = [tuple([-val + (phi * (val >= (phi >> 1))) for val in lig]) for lig in M1_or_B1]

		montgomery_convert_to_mns = montgomery_convert_to_mns_base
		amns_montg_mult = amns_montg_mult_base
		rho_div_convert_to_mns = rho_div_convert_to_mns_base

		# We cut up B and B1
		Bhi = [[elem >> 64 for elem in lig] for lig in B]
		Blo = [[elem % (2**64) for elem in lig] for lig in B]
		B1hi = [[elem >> 64 for elem in lig] for lig in B1]
		B1lo = [[elem % (2**64) for elem in lig] for lig in B1]

		print(f"static const uint64_t Blo[N][N] = {{\n{str(Blo)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', 'u},')}\n\t\t}},")
		print(f"\tB1lo[N][N] = {{\n{str(B1lo)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', 'u},')}\n\t\t}};")

		print(f"static const int64_t Bhi[N][N] = {{\n{str(Bhi)[1:-1].replace('[',  '{').replace(']', '}')}\n\t\t}},")
		print(f"\tB1hi[N][N] = {{\n{str(B1hi)[1:-1].replace('[',  '{').replace(']', '}')}\n\t\t}};")

		print("static const __int128 B1[N][N] = {")
		for i in range(n):
			print("\t{ ")
			for j in range(n - 1):
				print(f"((__int128) LOW({B1hi[i][j]}) << 64) | {B1lo[i][j]}u, ", end="")
			print(f"((__int128) LOW({B1hi[i][-1]}) << 64) | {B1lo[i][-1]}u }}", end="")
			if i != n - 1:
				print(",")
			else:
				print("\n\t\t};")

		# From here on, precalc functions
		print("""
static inline void m1_or_b1_mns128_mod_mult_ext_red_pre(unsigned __int128* restrict Rlo,
	const restrict poly128 A)
{""")

		for i in range(n):
			print(f"Rlo[{i}] = (__int128)", end="")
			for j in range(n):
				print(f" ((__int128) A->lo[{j}] * {B1lo[j][i]}u) + ((__int128) " +
					f"(LOW(A->lo[{j}] * {B1hi[j][i]}) + LOW(A->hi[{j}] * " +
					f"{B1lo[j][i]}u)) << 64)", end="")
				if j != n - 1:
					print(" +", end= "")

			print(";")

		print("}\n")

		print("""
static inline void m_or_b_mns128_mod_mult_ext_red_pre(__int128* restrict Rhi, 
	unsigned __int128* restrict Rlo, const restrict poly128 A)
{
	unsigned __int128 A0B0, A1B0, A0B1, tmplo;
	__int128 A1B1, aux2, aux3;
""")

		Bt = [Blo, Bhi]

		for i in range(n):
			for j in range(n):
				for k in range(3, -1, -1):
					print(f"\tA{k//2}B{k&1} = {At[k//2]}{j}]) * {Bt[k&1][j][i]}" +
						f"{'' if k&1 else 'u'};")
				print(auxstring)
				print(f"""\tRhi[{i}] += (__int128) aux2 + A1B1 +
		add_overflow(Rlo + {i}, tmplo);\n""")

		print("}\n")

	else:
		raise ValueError("Invalid parameter for M or B")

	# This one for the multiplication of A by B.
	print("""
static inline void mns128_mod_mult_ext_red_pre(__int128* restrict Rhi,
	unsigned __int128* restrict Rlo, const restrict poly128 A,
	const restrict poly128 B)
{
	unsigned __int128 A0B0, A1B0_A0B1, tmplo;
	__int128 A1B1;
""")

	At = ["(__int128) A->lo[", "(__int128) A->hi["]
	Bt = ["B->lo[", "B->hi["]

	for i in range(n):
		for j in range(1, n - i):
			print(f"\tA1B1 = {At[1]}{i + j}] * {Bt[1]}{n - j}] * LAMBDA;")
			print(f"\tA0B0 = {At[0]}{i + j}] * {Bt[0]}{n - j}] * LAMBDA;")
			print(f"\tA1B0_A0B1 = (__int128) ((__int128) ({At[1]}{i + j}] *")
			print(f"\t\t{Bt[0]}{n - j}] * LAMBDA) + ({At[0]}{i + j}] *")
			print(f"\t\t{Bt[1]}{n - j}] * LAMBDA)) + HIGH(A0B0);")
			print(f"tmplo = (__int128) LOW(A0B0) | ((__int128)A1B0_A0B1 << 64);")
			print(f"""\tRhi[{i}] += (__int128) A1B1 + HIGH(A1B0_A0B1) +
		add_overflow(Rlo + {i}, tmplo);\n""")

		for j in range(0, i + 1):
			print(f"\tA1B1 = {At[1]}{j}] * {Bt[1]}{i - j}];")
			print(f"\tA0B0 = {At[0]}{j}] * {Bt[0]}{i - j}];")
			print(f"\tA1B0_A0B1 = (__int128) ((__int128) ({At[1]}{j}] *")
			print(f"\t\t{Bt[0]}{i - j}]) + ({At[0]}{j}] *")
			print(f"\t\t{Bt[1]}{i - j}])) + HIGH(A0B0);")
			print(f"tmplo = (__int128) LOW(A0B0) | ((__int128)A1B0_A0B1 << 64);")
			print(f"""\tRhi[{i}] += (__int128) A1B1 + HIGH(A1B0_A0B1) +
		add_overflow(Rlo + {i}, tmplo);\n""")

	print("}\n")

	# We then get the Pi
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [0] * n
	Prho = montgomery_convert_to_mns(rho * phi, p, n, gamma, rho, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
	Pstk = montgomery_convert_to_mns(phi**2, p, n, gamma, rho, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
	for i in range(n):
		Pi[i] = Pstk
		Pstk = amns_montg_mult(Pstk, Prho, p, n, gamma, rho, lam, phi, M_or_B, M1_or_B1)
	print("\nstatic const uint64_t __Pilo__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str([int(elem) % 2**64 for elem in Pi[i]])[1:-1].replace(",", "u,") + "u},")
	print("\t\t{" + str([int(elem) % 2**64 for elem in Pi[-1]])[1:-1].replace(",", "u,") + "u}\n\t};\n")
	print("\nstatic const int64_t __Pihi__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str([int(elem) >> 64 for elem in Pi[i]])[1:-1] + "},")
	print("\t\t{" + str([int(elem) >> 64 for elem in Pi[-1]])[1:-1] + "}\n\t};\n")

	# We now transcribe the value of P
	tmp = convert_to_int_tabs(p)
	print("static _poly __P__ = { .deg = " + str(len(tmp)) + ",")
	tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
	print("\t\t.t = (int64_t[]) {" + tmp + "} },")

	# Powers of gamma next
	print("\tGi[] = { ", end="")
	g = gamma
	for i in range(1, n):
		tmp = convert_to_int_tabs(int(g))
		if i != 1:
			print("\t", end="")
		print("{ .deg = " + str(len(tmp)) + ",")
		print("\t\t.t = (int64_t[]) {", end="")
		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
		print(tmp + "} }", end="")
		if i != n - 1:
			print(",")
		g = g * gamma % p
	print("};")

	theta = rho_div_convert_to_mns(1, p, n, gamma, rho, lam, phi, M_or_B, M1_or_B1, Pi)
	tmphi = str([hex(int(elem) >> 64) for elem in theta])[1:-1].replace("'", "")
	tmplo = str([hex(int(elem) % 2**64) for elem in theta])[1:-1].replace("'", "")
	print(f"\nstatic _poly128 __theta__ = {{ .deg = {n},")
	print(f"\t.hi = (int64_t[]) {{ {tmphi} }},")
	print(f"\t.lo = (uint64_t[]) {{ {tmplo} }} }};")

	print("#endif")


if __name__ == "__main__":
#	p = 135235643069960614055763147653064061503447506836195743621526670176368184234720875133770858820476536434488752158042109408722131110950494765999584602783171012937442890496568352400298953673566640901275488915484314041362810104082501296367785554838617047203051373164870703338765011558032443485283660293658396373031
#	n = 9
#	gamma = 31655366034728588078624438514615301544129522139330836045091327852751402305472801538659044715201150659916554681038263189819329423658063171106562974285532277781484292697032498794895427121795477802587356170442668003146237461337340801098373462588444085927763403266750253334246732747899661111138239398839213291694
#	lam = 2
#	phi = 2**128
	if len(sys.argv) == 1:
		do_precalcs(p, n, gamma, lam, rho, M_or_B, M1_or_B1)
	else:
		try:
			psize = int(sys.argv[1])
			if int(str(psize) + "128") not in pmnsdicts.keys():
				print("Prime size not handled")
				exit()
			if len(sys.argv) > 2:
				index = int(sys.argv[2])
			else:
				index = 0
			index = primesdict[psize][index]
			pmns = pmnsdicts[int(str(psize) + "128")][index]
			do_precalcs(*pmns)
		except ValueError:
			print("Invalid syntax: [psize] [index]")
		except KeyError:
			print("Index value invalid")
	


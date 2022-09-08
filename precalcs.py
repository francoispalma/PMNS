import sys
from math import ceil

from ops import list_to_poly, montgomery_like_coefficient_reduction, horner_modulo, amns_montg_mult
from convert import montgomery_convert_to_mns, rho_div_convert_to_mns
from pyparams import p, n, gamma, lam, phi, rho, M, M1
from commonpmns import pmnsdicts, primesdict

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, rho, M, M1):
	print("#ifndef PMNS_PARAMS_H_INCLUDED\n#define PMNS_PARAMS_H_INCLUDED\n")

	print("#define RHO", rho)
	rho = 2**rho

	print("#define N", str(n) + "\n#define LAMBDA", str(lam) +"\n")

	# We print
	print("static const int64_t M[] = {" + str(M)[1:-1] + "},")
	print("\tM1[] = {" + str(M1)[1:-1] + "},")
	ML = [elem * lam for elem in M]
	print("\tMLambda[] = {" + str(ML)[1:-1] + "},")
	M1L = [(elem * lam) % phi for elem in M1]
	M1L = [M1L[i] - phi if M1L[i] >= (phi >> 1) else M1L[i] for i in range(n)]
	print("\tM1Lambda[] = {" + str(M1L)[1:-1] + "},")

	# We then get the Pi and print it
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), p, n, gamma, rho, lam, phi, M, M1, phinmoinsun) for i in range(n)]
	print("\t__Pi__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str(Pi[i])[1:-1] + "},")
	print("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n")

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

	theta = rho_div_convert_to_mns(1, p, n, gamma, rho, lam, phi, M, M1, Pi)
	tmp = str([hex(elem) for elem in theta])[1:-1].replace("'", "")
	print(f"\nstatic _poly __theta__ = {{ .deg = {n},")
	print(f"\t.t = (int64_t[]) {{ {tmp} }} }};")

	# From here on, precalc functions, first with the multiplication by M1
	print("""
static inline void m1_mns_mod_mult_ext_red_pre(int64_t* restrict R,
	const restrict poly A)
{
""")

	for i in range(n):
		print(f"R[{i}] = (uint64_t)", end="")
		for j in range(1, n - i):
			print(f" ((uint64_t) A->t[{i + j}] * {M1L[n - j]})", end="")
			print(" +", end= "")

		for j in range(0, i + 1):
			print(f" ((uint64_t) A->t[{j}] * {M1[i - j]})", end="")
			if j != i:
				print(" +", end= "")

		print(";")

	print("}\n")

	# Next is the multiplication by M
	print("""
static inline void m_mns_mod_mult_ext_red_pre(__int128* restrict R,
	const restrict poly A)
{
""")

	for i in range(n):
		print(f"R[{i}] = (__int128)", end="")
		for j in range(1, n - i):
			print(f" ((__int128) ((uint64_t) A->t[{i + j}]) * {ML[n - j]})", end="")
			print(" +", end= "")

		for j in range(0, i + 1):
			print(f" ((__int128) ((uint64_t) A->t[{j}]) * {M[i - j]})", end="")
			if j != i:
				print(" +", end= "")

		print(";")

	print("}\n")

	# This one for the multiplication of A by B.
	print("""
static inline void mns_mod_mult_ext_red_pre(__int128* restrict R,
	const restrict poly A, const restrict poly B)
{
""")

	for i in range(n):
		print(f"R[{i}] = (__int128)", end="")
		for j in range(1, n - i):
			print(f" ((__int128) A->t[{i + j}] * B->t[{n - j}] * LAMBDA)", end="")
			print(" +", end= "")

		for j in range(0, i + 1):
			print(f" ((__int128) A->t[{j}] * B->t[{i - j}])", end="")
			if j != i:
				print(" +", end= "")

		print(";")

	print("}\n")

	print("#endif")

if __name__ == "__main__":
#	p = 5524969863260095610495186344939419738027465949850580467176395575832917506871951737255621939342449907372936940924410124929668828406973361712220148691590192943
#	n = 12
#	gamma = 4636652768285458569062878611317602841499088282530798143147640460012527948722544350531006403703239417663454147165839392083945105225074778235517366830265723868
#	lam = -3
	if len(sys.argv) == 1:
		do_precalcs(p, n, gamma, lam, rho, M, M1)
	else:
		try:
			psize = int(sys.argv[1])
			if psize not in pmnsdicts.keys():
				print("Prime size not handled")
				exit()
			if len(sys.argv) > 2:
				index = int(sys.argv[2])
			else:
				index = 0
			index = primesdict[psize][index]
			pmns = pmnsdicts[psize][index]
			do_precalcs(*pmns)
		except ValueError:
			print("Invalid syntax: [psize] [index]")
		except KeyError:
			print("Index value invalid")

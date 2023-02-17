import sys
from math import ceil
from sage.all import matrix

from ops import horner_modulo, amns_montg_mult as amns_montg_mult_poly, amns_montg_mult_base
from convert import montgomery_convert_to_mns as montgomery_convert_to_mns_poly, montgomery_convert_to_mns_base, rho_div_convert_to_mns as rho_div_convert_to_mns_poly, rho_div_convert_to_mns_base
from pyparams import p, n, gamma, lam, phi, rho, M_or_B, M1_or_B1
from commonpmns import pmnsdicts, primesdict

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, rho, M_or_B, M1_or_B1):
	print("#ifndef PMNS_PARAMS_H_INCLUDED\n#define PMNS_PARAMS_H_INCLUDED\n")

	print("#define RHO", rho)
	rho = 2**rho

	print(f"#define N {n}")
	print('#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll ', end="")
	if n > 50:
		print('16")')
	else:
		print(f'{n}")')
	if type(lam) == int:
		print(f"#define LAMBDA {lam}")
	else:
		print(f"#define LENEXTPOLY {len(lam)}")
		print(f"static const int8_t EXTPOLY[{len(lam)}] = {{ {str(lam)[1:-1]} }};")

	# We determine if we're using the base matrix or a polynomial
	if type(M_or_B[0]) == int:
		print("#define M_or_B_is_M\n")
		M = M_or_B
		M = M + [0] * (n - len(M))
		M_or_B = M
		M1 = M1_or_B1

		montgomery_convert_to_mns = montgomery_convert_to_mns_poly
		amns_montg_mult = amns_montg_mult_poly
		rho_div_convert_to_mns = rho_div_convert_to_mns_poly

		# We print
		print("static const int64_t M[N] = {" + str(M)[1:-1] + "},")
		print("\tM1[N] = {" + str(M1)[1:-1] + "};")
#		ML = [elem * lam for elem in M]
#		print("\tMLambda[N] = {" + str(ML)[1:-1] + "},")
#		M1L = [(elem * lam) % phi for elem in M1]
#		M1L = [M1L[i] - phi if M1L[i] >= (phi >> 1) else M1L[i] for i in range(n)]
#		print("\tM1Lambda[N] = {" + str(M1L)[1:-1] + "};")

	elif type(M_or_B[0]) == list or type(M_or_B[0]) == tuple:
		print("#define M_or_B_is_B\n")
		B = M_or_B
		B1 = [tuple([-val + (phi * (val >= (phi >> 1))) for val in lig]) for lig in M1_or_B1]

		print(f"static const int64_t B[N][N] = {{\n{str(B)[1:-1].replace('(',  '{').replace(')', '}').replace('[',  '{').replace(']', '}')}\n\t\t}},")
		print(f"\tB1[N][N] = {{\n{str(B1)[1:-1].replace('(',  '{').replace(')', '}').replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}};")

		amns_montg_mult = amns_montg_mult_base
		montgomery_convert_to_mns = montgomery_convert_to_mns_base
		rho_div_convert_to_mns = rho_div_convert_to_mns_base

	# We then get the Pi and print it
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [0] * n
	Prho = montgomery_convert_to_mns(rho * phi, p, n, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
	Pstk = montgomery_convert_to_mns(phi**2, p, n, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
	for i in range(n):
		Pi[i] = Pstk
		Pstk = amns_montg_mult(Pstk, Prho, n, lam, phi, M_or_B, M1_or_B1)
	print("\nstatic const int64_t __Pi__[N][N] = {")
	for i in range(len(Pi) - 1):
		print("\t\t{" + str(Pi[i])[1:-1] + "},")
	print("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n")

	theta = rho_div_convert_to_mns(1, n, rho, lam, phi, M_or_B, M1_or_B1, Pi)
	tmp = str([hex(elem) for elem in theta])[1:-1].replace("'", "")
	print(f"\nstatic _poly __theta__ = {{ .deg = {n},")
	print(f"\t.t = (int64_t[]) {{ {tmp} }} }};")

	# We now transcribe the value of P
	tmp = convert_to_int_tabs(p)
	print("static _mpnum __P__ = { .deg = " + str(len(tmp)) + ",")
	print("\t\t.sign = 1,")
	tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
	print("\t\t.t = (uint64_t[]) {" + tmp + "} },")

	# Powers of gamma next
	print("\tGi[] = { ", end="")
	g = gamma
	for i in range(1, n):
		tmp = convert_to_int_tabs(int(g))
		if i != 1:
			print("\t", end="")
		print("{ .deg = " + str(len(tmp)) + ",")
		print("\t\t.sign = 1,")
		print("\t\t.t = (uint64_t[]) {", end="")
		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
		print(tmp + "} }", end="")
		if i != n - 1:
			print(",")
		g = g * gamma % p
	print("};")

	print("#endif")

if __name__ == "__main__":
	if len(sys.argv) == 1:
		do_precalcs(p, n, gamma, lam, rho, M_or_B, M1_or_B1)
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

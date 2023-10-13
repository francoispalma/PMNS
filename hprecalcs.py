import sys
from math import ceil, log2
from sage.all import matrix
from numpy import count_nonzero

from ops import horner_modulo, amns_montg_mult_base
from convert import montgomery_convert_to_mns_base, rho_div_convert_to_mns_base
from findm import findm
from generated.hpmns25664n5 import pmnsdict
#from generated.pmns51264n9 import pmnsdict
#from generated.pmns52164n9 import pmnsdict
#from generated.pmns41464n7 import pmnsdict
primes = list(pmnsdict.keys())
p, n, gamma, lam, rho, B, B1 = pmnsdict[primes[0]]
phi = 2**64

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, rho, B, B1):
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
		if abs(B[0][0]) != 1:
			M, M1 = findm(p, n, gamma, lam, rho, B, B1, phi)
			print(f"#define GAMMA {gamma}")
			if count_nonzero(M1) != 2:
				#lastcol = [(pow(p, -1, phi)*gamma**i) % phi for i in range(n)]
				print(f"static const uint64_t lastcol[{n}] = {{", end="")
				for i in range(n):
				 print(str((pow(p, -1, phi)*gamma**i) % phi) + "u, ")
				print("};\n")
			else:
				print("#define HOLLOWM1")
				print(f"#define GAMMALAMM1 {M1[-2]}")
				print(f"#define ONELAMM1 {M1[-1]}")
		else:
			M, M1 = B[0], B1[0]
			print("#define TWOTXMONE\n")
			print(f"#define TWOT {M[1]}")
			print(f"#define TWOTLAM {M[1] * lam}")
	else:
		print(f"#define BINOMIAL_A {lam['a']}")
		print(f"#define BINOMIAL_B {-lam['b']}")
		
		print(f"\n#define BINOMIAL_TWOT {B[1][1]}")
		print(f"#define BINOMIAL_TWOTOVERTWO {-B[0][0]}")


#	# We then get the Pi and print it
#	phinmoinsun = pow(phi, n - 1, p)
#	Pi = [0] * n
#	Prho = montgomery_convert_to_mns(rho * phi, p, n, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
#	Pstk = montgomery_convert_to_mns(phi**2, p, n, lam, phi, M_or_B, M1_or_B1, phinmoinsun)
#	for i in range(n):
#		Pi[i] = Pstk
#		Pstk = amns_montg_mult(Pstk, Prho, n, lam, phi, M_or_B, M1_or_B1)
#	print("\nstatic const int64_t __Pi__[N][N] = {")
#	for i in range(len(Pi) - 1):
#		print("\t\t{" + str(Pi[i])[1:-1] + "},")
#	print("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n")

#	theta = rho_div_convert_to_mns(1, n, rho, lam, phi, M_or_B, M1_or_B1, Pi)
#	tmp = str([hex(elem) for elem in theta])[1:-1].replace("'", "")
#	print(f"\nstatic _poly __theta__ = {{ .deg = {n},")
#	print(f"\t.t = (int64_t[]) {{ {tmp} }} }};")

#	# We now transcribe the value of P
#	tmp = convert_to_int_tabs(p)
#	print("static _mpnum __P__ = { .deg = " + str(len(tmp)) + ",")
#	print("\t\t.sign = 1,")
#	tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
#	print("\t\t.t = (uint64_t[]) {" + tmp + "} },")

#	# Powers of gamma next
#	print("\tGi[] = { ", end="")
#	g = gamma
#	for i in range(1, n):
#		tmp = convert_to_int_tabs(int(g))
#		if i != 1:
#			print("\t", end="")
#		print("{ .deg = " + str(len(tmp)) + ",")
#		print("\t\t.sign = 1,")
#		print("\t\t.t = (uint64_t[]) {", end="")
#		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
#		print(tmp + "} }", end="")
#		if i != n - 1:
#			print(",")
#		g = g * gamma % p
#	print("};")

	print(f"\n#define GMPLIMB {ceil(log2(p) - 1)//64 + 1}")
	print("\n#endif")

if __name__ == "__main__":
	do_precalcs(p, n, gamma, lam, rho, B, B1)

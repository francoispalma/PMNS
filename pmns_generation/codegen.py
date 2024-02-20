#!/bin/python3

import sys
from ast import literal_eval

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def exact_conversion_to_pmns(inp, p, n, phi, G, G1, phin):
	vecmat = lambda vec, mat: [sum([vec[j] * mat[j][i] for j in range(len(vec))]) for i in range(len(mat))]
	modphi = lambda vec: [(elem % phi) - phi * ((elem%phi) > phi/2) for elem in vec]
	plusvec = lambda vec1, vec2: [vec1[i] + vec2[i] for i in range(len(vec1))]
	slashphi = lambda vec: [elem//phi for elem in vec]
	Gmont_like = lambda vec: slashphi(plusvec(vec, vecmat(modphi(vecmat(vec, G1)), G)))
	A = [inp * phin % p] + [0] * (n-1)
	for i in range(n):
		A = Gmont_like(A)
	return A

def do_precalcs(p, n, gamma, E, rho, G, G1, phi, delta):
	with open("params.h", "w+") as out:
		out.write("#ifndef PMNS_PARAMS_H_INCLUDED\n#define PMNS_PARAMS_H_INCLUDED\n\n")
		out.write("#define RHO " + str(rho) + "\n")
		out.write(f"#define N {n}\n")
		if type(E) == int:
			out.write(f"#define LAMBDA {E}\n")
		elif type(E) == list:
			out.write(f"#define LENEXTPOLY {len(E)}\n")
			out.write(f"static const int8_t EXTPOLY[{len(E)}] = {{ {str(E)[1:-1]} }};\n")
		else:
			out.write(f"#define BINOMIAL_A {E['a']}\n")
			out.write(f"#define LAMBDA {-E['b']}\n")
		if phi != 64:
			out.write(f"#define PHI {phi}\n")
		phi = 1<<phi

		out.write("\n")

		out.write(f"static const int64_t G[N][N] = {{\n{str(G)[1:-1].replace('(',  '{').replace(')', '}').replace('[',  '{').replace(']', '}')}\n\t\t}};\n")
		out.write(f"static const int64_t G1[N][N] = {{\n{str(G1)[1:-1].replace('(',  '{').replace(')', '}').replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}};\n")

		# We then get the Pi and print it
		phin = pow(phi, n, p)
		Theta = rho.bit_length()
		while (1<<(Theta*n) < p):
			Theta += 1
		out.write("#define THETA " + str(Theta) + "\n")
		
		Pi = [0] * n
		Tau = phi**2 % p
		if type(E) == dict:
			Tau = Tau * int(pow(E['a'], -1, p)) % p
		for i in range(n-1):
			Pi[i] = exact_conversion_to_pmns(Tau, p, n, phi, G, G1, phin)
			Tau = Tau * 2**Theta % p
		Pi[n-1] = exact_conversion_to_pmns(Tau, p, n, phi, G, G1, phin)
		out.write("\nstatic const int64_t __Pi__[N][N] = {\n")
		for i in range(len(Pi) - 1):
			out.write("\t\t{" + str(Pi[i])[1:-1] + "},\n")
		out.write("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n\n")

		theta = exact_conversion_to_pmns(phi, p, n, phi, G, G1, phin)
		tmp = str([hex(elem) for elem in theta])[1:-1].replace("'", "\n")
		out.write(f"\nstatic _poly __theta__ = {{ .deg = {n},\n")
		out.write(f"\t.t = (int64_t[]) {{ {tmp} }} }};\n")

		# We now transcribe the value of P
		tmp = convert_to_int_tabs(p)
		out.write("static _mpnum __P__ = { .deg = " + str(len(tmp)) + ",\n")
		out.write("\t\t.sign = 1,\n")
		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "\n")
		out.write("\t\t.t = (uint64_t[]) {" + tmp + "} },\n")

		# Powers of gamma next
#		out.write("\tGi[] = { ", end="")
#		g = gamma
#		for i in range(1, n):
#			tmp = convert_to_int_tabs(int(g))
#			if i != 1:
#				out.write("\t", end="")
#			out.write("{ .deg = " + str(len(tmp)) + ",\n")
#			out.write("\t\t.sign = 1,\n")
#			out.write("\t\t.t = (uint64_t[]) {", end="")
#			tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "")
#			out.write(tmp + "} }", end="")
#			if i != n - 1:
#				out.write(",\n")
#			g = g * gamma % p
		tmp = convert_to_int_tabs(gamma)
		out.write("\tGi[] = { { .deg = " + str(len(tmp)) + ",\n")
		out.write("\t\t.sign = 1,\n")
		tmp = str([hex(elem) for elem in tmp])[1:-1].replace("'", "\n")
		out.write("\t\t.t = (uint64_t[]) {" + tmp + "} } };\n")
		out.write("#endif\n")

if __name__ == "__main__":
	if len(sys.argv) >= 2:
		try:
			with open(sys.argv[1]) as inpf:
				try:
					p, n, gamma, E, rho, G, G1 = literal_eval(inpf.readlines()[0])
					do_precalcs(p, n, gamma, E, rho, G, G1, phi=64, delta=0)
				except Exception as e:
					print("Invalid file format", e)
		except FileNotFoundError as e:
			print("File not found:", sys.argv[1])
		except Exception as e:
			print("Error:", e)
	else:
		print("oops")


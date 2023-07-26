from sage.all import matrix

from generated.pmns8192256 import pmnsdict

p, n, gamma, lam, rho, B, B1 = pmnsdict[list(pmnsdict.keys())[0]]

phi = 2**256

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def do_precalcs(p, n, gamma, lam, rho, M_or_B, M1_or_B1):
	print("#ifndef PMNS_PARAMS256_H_INCLUDED\n#define PMNS_PARAMS256_H_INCLUDED\n")

	print("#define RHO", rho)
	rho = 2**rho

	print(f"#define N {n}")
	print(f"#define LAMBDA {lam}\n")
	
	B = M_or_B
	B1 = [tuple([-val + (phi * (val >= (phi >> 1))) for val in lig]) for lig in M1_or_B1]

	# We cut up B and B1
	Bhi = [[elem >> 192 for elem in lig] for lig in B]
	Bmidhi = [[(elem >> 128) % (2**64) for elem in lig] for lig in B]
	Bmidlo = [[(elem >> 64) % (2**64) for elem in lig] for lig in B]
	Blo = [[elem % (2**64) for elem in lig] for lig in B]
	B1lo = [[elem % (2**64) for elem in lig] for lig in B1]
	
	
	print(f"static const uint64_t Blo[N][N] = {{\n{str(Blo)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}},")
	print(f"Bmidlo[N][N] = {{\n{str(Bmidlo)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}},")
	print(f"Bmidhi[N][N] = {{\n{str(Bmidhi)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}},")
	print(f"\tB1lo[N][N] = {{\n{str(B1lo)[1:-1].replace('[',  '{').replace(']', '}').replace(', ', 'u, ').replace('}u,', '},').replace('}', 'u}')}\n\t\t}};")

	print(f"static const int64_t Bhi[N][N] = {{\n{str(Bhi)[1:-1].replace('[',  '{').replace(']', '}')}\n\t\t}};")

	print("#endif")


if __name__ == "__main__":
	do_precalcs(p, n, gamma, lam, rho, B, B1)


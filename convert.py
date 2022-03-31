from random import randrange
from time import process_time

from ops import montgomery_like_coefficient_reduction, horner_modulo, naive_convert_to_mns, amns_montg_mult
from proof import amns, p, n, gamma, rho, phi, M, M1

def norme_infinie(L):
	return max([abs(elem) for elem in L])

def naive_convert_to_mns_and_red(a, p, n, gamma, rho, lam, phi, M, M1, phisquared):
	A = naive_convert_to_mns(a, p, n, gamma, rho, lam)
	U = amns_montg_mult(A, phisquared, p, n, gamma, rho, lam, phi, M, M1)
	return U

def montgomery_convert_to_mns(a, p, n, gamma, rho, lam, phi, M, M1, tau):
	alpha = (a * tau) % p
	A = [0] * n
	A[0] = alpha
	for _ in range(n - 1):
		A = montgomery_like_coefficient_reduction(A, p, n, gamma, rho, lam, phi, M, M1)
	return A

def rho_div_convert_to_mns(a, p, n, gamma, rho, lam, phi, M, M1, Pi):
	t = convert_to_rho_base(a, n, rho)
#	print("rho div tab")
#	print([hex(elem) for elem in t])
#	print()
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
#	print("after Pi")
#	print([hex(elem) if elem >= 0 else hex(2**256 + elem) for elem in U])
#	print()
#	print("HAAAAAAAAAA")
#	print("M:")
#	print([hex(elem) for elem in M])
#	print("M1:")
#	print([hex(elem) for elem in M1])
#	print("Pi:")
#	print([[hex(Pi[i][j]) for j in range(n)] for i in range(n)])
	A = montgomery_like_coefficient_reduction(U, p, n, gamma, rho, lam, phi, M, M1)
	return A

def convert_to_rho_base(a, n, rho):
	t = [0] * n
	a1 = a
	RHO = rho.bit_length() - 1
	for i in range(n):
		t[i] = a1 & (rho - 1)
		a1 >>= RHO
	return t

if __name__ == "__main__":
	phisquared = [0, 0, 0, 512, 0]
	phinmoinsun = pow(phi, n - 1, p)
	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), *amns, phi, M, M1, phinmoinsun) for i in range(n)]
	sum1 = 0
	sum2 = 0
	sum3 = 0
#	sum4 = 0
	for _ in range(10000):
		A = randrange(p)
		c = process_time()
		A1 = naive_convert_to_mns_and_red(A, *amns, phi, M, M1, phisquared)
		sum1 += process_time() - c
		c = process_time()
		A2 = montgomery_convert_to_mns(A, *amns, phi, M, M1, phinmoinsun)
		sum2 += process_time() - c
		c = process_time()
		A3 = rho_div_convert_to_mns(A, *amns, phi, M, M1, Pi)
		sum3 += process_time() - c
#		c = process_time()
#		A4 = naive_convert_to_mns(A * phi, *amns)
#		sum4 += process_time() - c
		if horner_modulo(A1, gamma, p) != (A * phi % p) or norme_infinie(A1) >= rho:
			print("Naive Error")
		if horner_modulo(A2, gamma, p) != A or norme_infinie(A2) >= rho:
			print("Montgom Error")
		if horner_modulo(A3, gamma, p) != (A * phi % p) or norme_infinie(A3) >= rho:
			print("Rhodiv Error")
	print(sum1)
	print(sum2)
	print(sum3)
#	print(sum4)
	
	# Convert from amns to binary
	a = 0x77f882926258fb5a293015e16fc961598939f9f328d4e316d02519d3f8d88412d787
	b = 0xb4399ccbab87f4f053d75a9dcc1c1fa8d2f4edd7bdf5eebc78fb4ea16a6fb02eb96d
	c = a * b
	A = rho_div_convert_to_mns(a, *amns, phi, M, M1, Pi)
	B = rho_div_convert_to_mns(b, *amns, phi, M, M1, Pi)
	C = amns_montg_mult(A, B, *amns, phi, M, M1)
	print([hex(elem) for elem in C])
	C = montgomery_like_coefficient_reduction(C, *amns, phi, M, M1)
	print([hex(elem) for elem in C])
	G1 = gamma**1 % p
	G2 = gamma**2 % p
	G3 = gamma**3 % p
	G4 = gamma**4 % p
	Gi = [1, G1, G2, G3, G4]
	cc = C[0]
	for i in range(1, n - 1):
		print(hex(cc))
		cc += C[i] * Gi[i]
	print(hex(cc))
	cc = cc % p
	print(hex(cc))
	
		

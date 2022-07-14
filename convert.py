from random import randrange
from time import process_time
from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice

from ops import montgomery_like_coefficient_reduction, horner_modulo, naive_convert_to_mns, amns_montg_mult, list_to_poly
from pyparams128 import p, n, gamma, rho, lam, M, M1, phi

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
	U = [0] * n
	for i in range(n):
		for j in range(n):
			U[j] += t[i] * Pi[i][j]
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
	pmns = (p, n, gamma, rho, lam, phi, M, M1)
	
	phinmoinsun = pow(phi, n - 1, p)
	
	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), *pmns, phinmoinsun) for i in range(n)]
	
	a = 0xeedad1f0fcbd8935abfccb2a2752c672db377b0575464c626bd4003dca25c7c80589ac8e48c02c2379c73d45d9558ba5e708e338ec0686a4ab1356e9ad34348859398e506895f6c321bc4f8d2e2f21783ff2e65030a968331c9829e17900a1af32f1b7a3131ab61c0ad9b9e98a23f932d2f8936bc32d06128b086e619a1e35c1
	b = 0xa3442cfd10dcb1d9a20a8541a8551a37d23bbe02dff3faeaf4fd5253975bfad4c759b9101bf4ae5197a18385f0658f5e436644638cd923ef490a7251402fdc0ef87dd61af1f0dbc9b81ef7de744b062a0d5df758eb1c1c41fd4cf01255aba467f1f673c4e5b20540f9c77ae0e430a323225346fdcc804912a09c6728b5837935
	
	#A = rho_div_convert_to_mns(a, *pmns, Pi)
	A = montgomery_convert_to_mns(a, *pmns, phinmoinsun)
	
	#B = rho_div_convert_to_mns(b, *pmns, Pi)
	B = montgomery_convert_to_mns(b, *pmns, phinmoinsun)
	
	print(A)
	print(hex(a))
	print(hex(horner_modulo(A, gamma, p) * 1 % p))
	print()
	print(B)
	print(hex(b))
	print(hex(horner_modulo(B, gamma, p) * 1 % p))

#if __name__ == "__main__":
#	phisquared = [0, 0, 0, 512, 0]
#	phinmoinsun = pow(phi, n - 1, p)
#	Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), *amns, phi, M, M1, phinmoinsun) for i in range(n)]
#	sum1 = 0
#	sum2 = 0
#	sum3 = 0
##	sum4 = 0
#	for _ in range(10000):
#		A = randrange(p)
#		c = process_time()
#		A1 = naive_convert_to_mns_and_red(A, *amns, phi, M, M1, phisquared)
#		sum1 += process_time() - c
#		c = process_time()
#		A2 = montgomery_convert_to_mns(A, *amns, phi, M, M1, phinmoinsun)
#		sum2 += process_time() - c
#		c = process_time()
#		A3 = rho_div_convert_to_mns(A, *amns, phi, M, M1, Pi)
#		sum3 += process_time() - c
##		c = process_time()
##		A4 = naive_convert_to_mns(A * phi, *amns)
##		sum4 += process_time() - c
#		if horner_modulo(A1, gamma, p) != (A * phi % p) or norme_infinie(A1) >= rho:
#			print("Naive Error")
#		if horner_modulo(A2, gamma, p) != A or norme_infinie(A2) >= rho:
#			print("Montgom Error")
#		if horner_modulo(A3, gamma, p) != (A * phi % p) or norme_infinie(A3) >= rho:
#			print("Rhodiv Error")
#	print(sum1)
#	print(sum2)
#	print(sum3)
##	print(sum4)
#	
#	# Convert from amns to binary
#	a = 0x77f882926258fb5a293015e16fc961598939f9f328d4e316d02519d3f8d88412d787
#	b = 0xb4399ccbab87f4f053d75a9dcc1c1fa8d2f4edd7bdf5eebc78fb4ea16a6fb02eb96d
#	c = a * b
#	A = rho_div_convert_to_mns(a, *amns, phi, M, M1, Pi)
#	B = rho_div_convert_to_mns(b, *amns, phi, M, M1, Pi)
#	C = amns_montg_mult(A, B, *amns, phi, M, M1)
#	print([hex(elem) for elem in C])
#	C = montgomery_like_coefficient_reduction(C, *amns, phi, M, M1)
#	print([hex(elem) for elem in C])
#	G1 = gamma**1 % p
#	G2 = gamma**2 % p
#	G3 = gamma**3 % p
#	G4 = gamma**4 % p
#	Gi = [1, G1, G2, G3, G4]
#	cc = C[0]
#	for i in range(1, n - 1):
#		print(hex(cc))
#		cc += C[i] * Gi[i]
#	print(hex(cc))
#	cc = cc % p
#	print(hex(cc))
	

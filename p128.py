from sage.all import matrix, ZZ, PolynomialRing
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

if __name__ == "__main__":
	p = 5524969863260095610495186344939419738027465949850580467176395575832917506871951737255621939342449907372936940924410124929668828406973361712220148691590192943
	n = 12
	gamma = 4636652768285458569062878611317602841499088282530798143147640460012527948722544350531006403703239417663454147165839392083945105225074778235517366830265723868
	#rho = 98
	lam = -3
	B = [[p if (i, j) == (0, 0) else pow(-gamma, i, p) if i != 0 and j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]
	B = list(IntegerLattice(matrix(ZZ, B)).LLL())
	RingPoly = PolynomialRing(ZZ, 'X')
	E = RingPoly("X^12 - 3")
	phi = 2**128
	w = 1 + (n - 1) * abs(lam)
	__tmp = int(2 * w * max(max(B)))
	rho = ceil(__tmp.bit_length())
	print("rho = 2**" + str(rho))
	rho = 2**rho
	for lig in B:
		for elem in lig:
			if abs(elem) > rho:
				print("Something wrong with LLL'd matrix")
		M = RingPoly(list_to_poly(lig))
		val, M1, soak = xgcd(M, E)
		val = ZZ(val)
		if val & 1:  # if val is even it won't have a gcd of 1 with phi.
			M1 = (M1 * ZZ(pow(val, -1, phi)) % phi)
			print(((M * M1) % E) % phi)
			print(M(gamma) % p == 0)
			break
	M = list(M)
	print("M")
	print(M)
	M1 = list(M1)
	print("M1")
	M1 = [(int(M1[i]) * -1) % phi for i in range(n)]
	M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]
	print(M1)

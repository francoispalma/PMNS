from sage.all import matrix, ZZ, PolynomialRing, xgcd
from sage.modules.free_module_integer import IntegerLattice
from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from pyparams import p, n, gamma, lam
from ops import list_to_poly, horner_modulo

RHO = 118
rho = 2**RHO
phi = 2**128

B = [[p if (i, j) == (0, 0) else -pow(gamma, i, p) if i != 0 and j == 0 else 1 if i == j else 0 for j in range(n)] for i in range(n)]

# We apply LLL to it
B = list(IntegerLattice(matrix(ZZ, B)).LLL())

# We define our external reduction polynomial
RingPoly = PolynomialRing(ZZ, 'X')
E = RingPoly("X^" + str(n) + " - (" + str(lam) + ")")

for lig in B:
	# We check if something went wrong
	for elem in lig:
		if abs(elem) > rho:
			print("Something wrong with LLL'd matrix")
			exit()

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
M1 = list(M1)

# We switch M1 from M^-1 to -M^-1
M1 = [(int(M1[i]) * -1) % phi for i in range(n)]

# We then reduce M1's coefficients
M1 = [int(M1[i]) - phi if M1[i] >= (phi >> 1) else M1[i] for i in range(n)]

phinmoinsun = pow(phi, n - 1, p)
Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), p, n, gamma, rho, lam, phi, M, M1, phinmoinsun) for i in range(n)]

a = 0x74ff560400d0105e6381e4f7cf22ba4a3d949bbe3b03e7ec1c8aebfb02a4dedf230eef099cd1ae78adf8f142cd70ed93122a5c48c5edcba658615fa2316994dce0c84e9e54c5ae9482acdc0ed6fae84eb7e83d94016d12452ad41369e33a53a676d539439488bdc8b3462c5579a432e8b579e8af9d5b2b0b8f37856fe2de7f30
print("HAAAAAAAAAA")
tmp = conv(a, p, n, gamma, rho, lam, phi, M, M1, Pi)

print(tmp)

print()

tmp = horner_modulo(tmp, gamma, p)
#print(tmp)

print()
print(tmp == a * phi % p)


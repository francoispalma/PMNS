from sage.all import ZZ, PolynomialRing, xgcd

RingPoly = PolynomialRing(ZZ, "X")

def findm(p, n, gamma, lam, rho, B, B1, PHI):
	E = RingPoly(f"X^{n}") - RingPoly(lam)
	# We then try to find a valid M
	for lig in B:
		# We set M as one line of the base matrix
		M = RingPoly(lig)
		# We then get the d and u of the au + bv = d from extended euclid
		val, M1, soak = xgcd(M, E)
		# We make sure to get out of any sage specific field
		try:
			val = ZZ(val)
		except TypeError:
			val = 0
		if val & 1:  # if val is even it won't have a gcd of 1 with phi
			M1 = (M1 * ZZ(pow(val, -1, PHI)) % PHI)
			# We check that M1 is properly constructed, if it is we don't look further
			if ((M * M1) % E) % PHI == 1 and M(gamma) % p == 0:
				break
	# We convert out of the polynomial Ring
	M = list(M)
	M1 = list(M1)
	# We switch M1 from M^-1 to -M^-1
	M1 = [(int(M1[i]) * -1) % PHI for i in range(n)]
	# We then reduce M1's coefficients
	M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else M1[i] for i in range(n)]
#	__tmp = int(2 * w * max([abs(elem) for elem in M]))
#	rho = __tmp.bit_length()
	return M, M1



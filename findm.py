from sage.all import ZZ, PolynomialRing, xgcd

RingPoly = PolynomialRing(ZZ, "X")

def infinite_norm_of_poly(poly):
	return max([abs(elem) for elem in poly])

def findm(p, n, gamma, lam, rho, B, B1, PHI):
	E = RingPoly(f"X^{n}") - RingPoly(lam)
	validMs = []
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
				validMs += [(list(M), list(M1))]
	if len(validMs) == 0:
		return []
	minMcouple = []
	minnorm = p
	for couple in validMs:
		polynorm = infinite_norm_of_poly(couple[0])
		if polynorm < minnorm:
			minnorm = polynorm
			minMcouple = couple
	# We convert out of the polynomial Ring
	M, M1 = minMcouple
	M = [int(elem) for elem in M]
	# We switch M1 from M^-1 to -M^-1
	M1 = [(int(M1[i]) * -1) % PHI for i in range(n)]
	# We then reduce M1's coefficients
	M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else M1[i] for i in range(n)]
	w = 2 * n - 1
	__tmp = int(2 * w * max([abs(elem) for elem in M]))
	rho = __tmp.bit_length()
	return M, M1



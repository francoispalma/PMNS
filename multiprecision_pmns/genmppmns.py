#!/usr/bin/env python3
import sys

from sage.all import ZZ, random_prime, GF, factor, matrix, round,xgcd,vector, ceil, gcd, exp, log

VERBOSE = True  # enables verbose mode
FAST = True  # if true skips finding polynomial roots through sage

s = 3  # parameter s
size = 8192  # prime size
n = 48  # parameter n
maxlambda = 2  # max value of lambda allowed in the generation
SIGMA = 64  # size of a memory register
PHI = s*SIGMA  # temporary phi for generation purposes
SANITY = 2**SIGMA  # sanity check at the end

def eprint(*args, **kwargs):
	if VERBOSE:
		print(*args, file=sys.stderr, **kwargs)

NBCHUNKS = s

if size//(SIGMA-ceil(log(size/(SIGMA*exp(1)),2))) > n * NBCHUNKS:
	eprint("Warning: n seems too small for this size")

def findminrow(mat, lam):
	"""Function that finds the row with the smallest 1-norm.
	"""
	norm1 = lambda V: sum([abs(elem) for elem in V])
	# The 1-norm of the companion matrix is computed with compnorm.
	compnorm = lambda X : norm1(vector(list(X)[1:]))*lam + abs(X[0])
	minM = []
	minnorm = abs(mat.det())
	mindex = -1
	i = 0
	E = ZZ["X"](f"X^{n} - {lam}")
	for lig in mat:
		tmp = list(lig)
		while True:
			M = ZZ["X"](tmp)
			val, M1, soak = xgcd(M, E)
			try:
				val = ZZ(val)
				if val & 1:  # if val is even it won't have a gcd of 1 with phi
					norm = compnorm(tmp)
					if norm < minnorm:
						minnorm = norm
						mindex = i
						minM = tmp
			except TypeError:
				pass
			if tmp[0]%lam:
				break
			tmp = tmp[1:] + [tmp[0]//lam]
		i += 1
	return mindex, minnorm, minM

if __name__ == "__main__":
	ZZX = ZZ["X"]
	chunks = NBCHUNKS
	cpt = 0
	notfound = True
	while notfound:
		eprint(f"\bCandidate number {cpt+1}", end="\r")
		p = random_prime(2**size,lbound=2**(size-1),proof=False)
		K = GF(p,proof=False)
		lam = 2
		# We try to find the nth root fast
		hdrb = gcd(p-1,n)
		if hdrb == 1:  # case GCD = 1
			unused, u, unused = xgcd(n, p-1)
			rts = [[int(pow(lam, u, p))]]
			if pow(rts[0][0], n, p) != lam:
				eprint("\nOops")
				eprint(p, n, lam, hdrk, rts[0][0])
				rts = []
		else:  # case GCD != 1
			while lam <= maxlambda and 2*(lam*(n-1) + 1)*round(n/exp(1))*ceil(2**((size-1)/n)) < 2**PHI:
				hdrk = 1
				while int(pow(lam,((p-1)//(hdrb**hdrk)),p)) == 1:
					hdrk += 1
				hdrk -= 1
				while hdrk > 0 and gcd((p-1)//(hdrb**hdrk), n) != 1:
					hdrk -= 1
				if hdrk > 0:
					hdru = pow(n, -1, (p-1)//(hdrb**hdrk))
					rts = [[int(pow(lam, int(hdru), p))]]
					if pow(rts[0][0], n, p) != lam:
						eprint("\nOops")
						eprint(p, n, lam, hdrk, rts[0][0])
						exit()
					break
				else:
					lam += 1
			else:  # If we do not find anything we ask sage
				rts = []
				lam = 2
				if not FAST:
					while lam <= maxlambda and 2*(lam*(n-1) + 1)*round(n/exp(1))*ceil(2**((size-1)/n)) < 2**PHI:
						E = ZZX(f"X^{n} - {lam}")
						rts = E.roots(ring=K)
						if rts:
							break
						lam += 1
		if not rts:  # We cannot find a PMNS with these parameters, we choose a new p
			cpt += 1
			continue
		# We have found a valid PMNS, now we construct it
		normin = p
		gammamin = -1
		for root in rts:  # We go through each gamma candidate
			gamma = root[0]
			G = matrix(ZZ,[[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]).LLL()
			n1G = findminrow(G, lam)[1]
			if n1G < normin:
				normin = n1G
				gammamin = gamma
		# We have found the gamma which gives the smallest 1-norm matrix
		w = (n-1) * abs(lam) + 1
		if 2*w*(normin - 1) < 2**PHI:
			rho = normin - 1
			phi = 2*w*rho
			philog2 = int(phi).bit_length()
			rhover2 = ceil(philog2/chunks)
			eprint("\nfound\n")
			if n*lam*rho**(ZZ(1)/NBCHUNKS) < SANITY:  # Last sanity check
				break
			else:
				eprint("bit big", round(n*lam*rho**(ZZ(1)/NBCHUNKS)/SANITY), rhover2, hdrb)
		else:
			eprint("\ntoo big:", int(2*w*(normin - 1)) >> PHI)
		cpt += 1
	# We have found a valid PMNS so we print
	print((p,n,gammamin,lam,normin-1))
	eprint(rhover2)


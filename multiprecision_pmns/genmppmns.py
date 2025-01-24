#!/usr/bin/env python3
import sys

from sage.all import ZZ, random_prime, GF, factor, matrix, round,xgcd,vector, ceil, gcd, exp, log

VERBOSE = True
FAST = True

NBCHUNKS = 4
PHI = NBCHUNKS*64
size = 8192
n = 144//4
maxlambda = 2

def eprint(*args, **kwargs):
	if VERBOSE:
		print(*args, file=sys.stderr, **kwargs)

if size//(64-ceil(log((size//64)/exp(1),2))) > n * NBCHUNKS:
	eprint("Warning: n seems too small for this size")

def findminrow(mat, lam):
	"""Function that finds the row with the smallest 1-norm.
	"""
	norm1 = lambda V: sum([abs(elem) for elem in V])
	# The 1-norm of the companion matrix is computed with compnorm.
	compnorm = lambda X : norm1(vector(list(X)[1:]))*lam + abs(X[0])
	minnorm = mat.det()
	mindex = -1
	i = 0
	for lig in mat:
		tmp = list(lig)
		# We can reduce the norm slightly by dividing by X mod E.
		while not(tmp[0]%lam):
			tmp = tmp[1:] + [tmp[0]//lam]
		norm = compnorm(tmp)
		if norm < minnorm:
			minnorm = norm
			mindex = i
		i += 1
	return mindex, minnorm

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
			phi = 2*w*(normin - 1)
			philog2 = int(phi).bit_length()
			rhover2 = ceil(philog2/chunks)
			eprint("\nfound\n")
			if 2*lam*2**rhover2 < 2**64:  # Last sanity check
				break
			else:
				eprint("bit big", (2*lam*2**rhover2)>>64, rhover2)
		else:
			eprint("\ntoo big:", int(2*w*(normin - 1)) >> PHI)
		cpt += 1
	# We have found a valid PMNS so we print
	print((p,n,gammamin,lam,normin-1))
	eprint(rhover2)


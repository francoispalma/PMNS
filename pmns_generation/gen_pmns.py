#!/usr/bin/python3

from sage.all import next_prime, next_probable_prime, previous_prime, GF, PolynomialRing, factor, matrix, ZZ, xgcd, random_prime, is_prime, vector, gcd, factor, pari, randrange
from math import log2
import sys
import datetime
import gc

Polydict = {
	"X^{n} - X^{n//2} + 1":(lambda n: 3*n/2, lambda n:[-1] + [0] * (n//2 - 1) + [1]), "X^{n} + X^{n//2} + 1":(lambda n: 3*n/2, lambda n:[-1] + [0] * (n//2 - 1) + [-1]),
	"X^{n} - X - 1":(lambda n: 2*n - 1, [1, 1]),"X^{n} - X + 1":(lambda n: 2*n - 1, [-1, 1]), "X^{n} + X - 1":(lambda n: 2*n - 1, [1, -1]), "X^{n} + X + 1":(lambda n: 2*n - 1, [-1, -1])
	}

VERBOSE = False
ZZX = ZZ["X"]

nthroot = lambda p, n: 2**(log2(p)/n)

def eprint(*args, **kwargs):
	if VERBOSE:
		print(*args, file=sys.stderr, **kwargs)


def findeasynthroot(prime, N):
	n = N
	rat = prime - 1
	while rat % n == 0:
		rat = rat//n
	while n % 2 == 0 and rat % 2 == 0 and n > 2:
		rat = rat // 2
		n = n // 2
	k = 1
	while k < n and (k * rat + 1) % n != 0:
		k += 1
	u = ((k * rat + 1)) // n
	if (pow(pow(2, u, prime), n, prime) == 2):
		return u
	else:
		return None

def get_phi_bound(n, w, nthrootofp, delta):
	return 2 * w * (1.02)**n * nthrootofp * ((1 + delta)**2) * (n**0.5)

def get_init_n(p, PHI, delta=0):
	init_n = p.bit_length()//int(log2(PHI)) + 1
	while PHI < get_phi_bound(init_n, 3*init_n/2, nthroot(p, init_n), delta):
		init_n += 1
	if PHI < get_phi_bound(init_n, 2*init_n - 1, nthroot(p, init_n), delta):
		init_n += 1
	return init_n


def checkpoly(p, rts, n, w, lam, nthrootofp, delta, PHI):
	E = ZZX(f"X^{n}") - ZZX(lam)
	if len(rts) == 0:
		return []

	for rt in rts:
		gamma = int(rt[0]) % p
		if gamma <= nthrootofp or gamma >= p - nthrootofp:
			continue
		if type(lam) == int and pow(gamma, n, p) != (lam % p):
			continue

		B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
		B = list(matrix(ZZ, B).LLL())
		n1B = max([sum([abs(B[i][j]) for i in range(len(B))]) for j in range(len(B))])
		two_double_u = round(2*w)
		if two_double_u*(n1B-1)*((1+delta)**2) <= PHI:
			B1 = list(matrix(B).inverse() % PHI)
			rho = n1B - 1
			return [p, n, gamma, lam, rho, B, B1]
	return []


def check_eval_strE(eval_strE, polK, p, n, w, lam, nthrootofp, delta, PHI):
	E = ZZX(eval_strE)
	highest_degree = n - len(factor(E))
	if PHI < get_phi_bound(n, w, nthroot(p, highest_degree), delta):
		return []
	try:
		rts = E.roots(ring=polK)
	except cypari2.handle_error.PariError:
		pari.allocatemem()
		rts = E.roots(ring=polK)
	check = checkpoly(p, rts, n, w, lam, nthrootofp, delta, PHI)
	if check:
		eprint(eval_strE)
	return check

def genpmns(p, PHI, init_n=None, amns_only=False, max_lambda=1<<64, delta=0):
	if init_n is None or init_n == 0:
		init_n = get_init_n(p, PHI, delta)
	K = None
	polK = None
	check = []
	n = max(init_n, 2)
	eprint("Initial n:", n)
	if amns_only:
		polydict = {}
	else:
		polydict = Polydict
	while True:
		nthrootofp = int(nthroot(p, n))
		#AMNS: X^n - lambda
		
		pgcd, u, v = xgcd(n, p - 1)
		d = pgcd
		if d == 1:
			check = checkpoly(p, [[int(pow(2, u, p))]], n, 2*n - 1, 2, nthrootofp, delta, PHI)
			if check:
				eprint(f"X^{n} - 2 found gcd 1")
				break
		u = findeasynthroot(p, n)
		if u:
			check = checkpoly(p, [[int(pow(2, u, p))]], n, 2*n - 1, 2, nthrootofp, delta, PHI)
			if check:
				eprint(f"X^{n} - 2 found easy")
				break
		lam = 2
		w = lambda n:lam * (n - 1) + 1
		if d != 1:
			while lam < max_lambda and PHI >= get_phi_bound(n, w(n), nthrootofp, delta):
				for lamb in [lam, -lam]:
					d = pgcd
					if pow(lamb, (p-1)//d, p) == 1:
						overbk = p - 1
						maxkgcd1 = p - 1
						while overbk % d == 0:
							overbk //= d
							maxkgcd1 = overbk if gcd(overbk, n) == pow(lamb, overbk, p) == 1 else maxkgcd1
						if maxkgcd1 != p - 1:
							soak, u, soak = xgcd(n, maxkgcd1)
							check = checkpoly(p, [[int(pow(lamb, u, p))]], n, w(n), lamb, nthrootofp, delta, PHI)
							if check:
								eprint(f"X^{n} - {lamb} found")
								break
						if check:
							break
				else:
					lam += 1
					continue
				break
		if check:
			break
		if K is None:
			eprint("Generating GF(p)", datetime.datetime.now().time())
			polK = GF(p, proof=False)
			eprint("factors", datetime.datetime.now().time())
		lam = 1
		while lam < max_lambda and PHI >= get_phi_bound(n, w(n), nthrootofp, delta):
			eprint("Checking lambda", lam, datetime.datetime.now().time())
			eval_strE = f"X^{n} - {lam}"
			check = check_eval_strE(eval_strE, polK, p, n, w(n), lam, nthrootofp, delta, PHI)
			if check:
				break
			eval_strE = f"X^{n} + {lam}"
			check = check_eval_strE(eval_strE, polK, p, n, w(n), -lam, nthrootofp, delta, PHI)
			if check:
				break
			lam += 1
		if check:
			break

		eprint("No amns")
		#PMNS
		for strE in polydict:
			if "//" in strE and n & 1:
				continue
			w, lam = polydict[strE]
			w = w(n)
			if callable(lam):
				lam = lam(n)
			eval_strE = eval(f'f"{strE}"')
			check = check_eval_strE(eval_strE, polK, p, n, w, lam, nthrootofp, delta, PHI)
			if check:
				#print(check)
				break
		else:
			n += 1
			eprint("Checking n =", n, "and freeing", gc.collect())
			if n > 2 * init_n:
				break
			continue
		break
	if check != []:
		return check
	eprint("finished", datetime.datetime.now().time())

def handle_opts(args, opts, optname, shortopt, longopt):
	if shortopt not in opts and longopt not in opts:
		return None
	if shortopt in opts:
		op = shortopt
	else:
		op = longopt
	try:
		retval = int(args[args.index(op) + 1])
	except IndexError:
		print("Error in use of option:", op)
		print("Not enough arguments")
		exit()
	except ValueError:
		print("Invalid value for", optname, args[args.index(op) + 1])
		exit()
	args.pop(args.index(op) + 1)
	return retval

if __name__ == "__main__":
	args = sys.argv.copy()
	if "-v" in args:
		VERBOSE = True
		args.pop(args.index("-v"))
	if "--verbose" in args:
		VERBOSE = True
		args.pop(args.index("--verbose"))
	amns_only = False
	if "--AMNS" in args:
		amns_only = True
		args.pop(args.index("--AMNS"))
	opts = [arg for arg in args if arg.startswith("-")]
	p = handle_opts(args, opts, "p", "-p", "--prime")
	if p is not None:
		try:
			assert p >= 2, f"Invalid value for p: {p}"
			assert is_prime(p), f"Value for p not prime: {p}"
		except Exception as e:
			print(e)
			exit()
	delta = handle_opts(args, opts, "delta", "-d", "--delta")
	if delta is not None:
		try:
			assert delta >= 0, f"Invalid value for delta: {delta}"
		except Exception as e:
			print(e)
			exit()
	else:
		delta = 0
	phi = handle_opts(args, opts, "phi", "-F", "--phi")
	if phi is not None:
		try:
			assert phi >= 0, f"Invalid value for phi: {phi}"
		except Exception as e:
			print(e)
			exit()
	else:
		phi = 64
	max_lambda = handle_opts(args, opts, "max_lambda", "-l", "--lambda")
	if max_lambda is not None:
		try:
			assert max_lambda >= 0, f"Invalid value for lambda: {max_lambda}"
		except Exception as e:
			print(e)
			exit()
	else:
		max_lambda = 1<<64
	arguments = [arg for arg in args if not arg.startswith("-")]
	if p is None and len(arguments) >= 2:
		try:
			psize = int(arguments[1])
		except ValueError:
			print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
			exit()
		if psize <= 1:
			print("Prime Size not handled")
			exit()
		if psize < 1664:
			p = int(random_prime(2**psize, lbound=2**(psize-1)))
		else:
			p = int(next_probable_prime(randrange(2**(psize-1), 2**psize)))
	if p:
		eprint()
		if amns_only:
			eprint("Generating AMNS", end=" ")
		else:
			eprint("Generating PMNS", end=" ")
		eprint("for p =", p)
		eprint("With PHI = 2^"+ str(phi))
		if delta:
			eprint("With delta =", delta)
		eprint()
		eprint("Output format [p, n, gamma, E, rho, B, B1]", end="")
		if amns_only:
			eprint()
		else:
			eprint(" with E as either an integer for AMNS or a list for non-AMNS")
		eprint("\n")
		soak, n, gamma, E, rho, B, B1 = genpmns(p, 2**phi, amns_only=amns_only, delta=delta, max_lambda=max_lambda)
		if type(E) == int:
			assert (gamma**n - E) % p == 0
		else:
			assert (gamma**n - ZZX(E)(gamma)) % p == 0
		print([p, n, gamma, E, rho, B, B1])
	else:
		print("Usage: ./gen_pmns.py [OPTION] [prime size]\n")
		print("\t--prime, -p: generate pmns for specific prime instead of a random one.")
		print("\t--phi, -F: sets the parameter PHI = 2^phi. Default value is phi=64.")
		print("\t--delta, -d: use specified delta for number of free additions.")
		print("\t--AMNS: turns generation into AMNS-only mode which is faster but will ignore other good PMNS shapes.")
		print("\t--lambda, -l: use specified value for maximum lambda allowed in AMNS.")
		print("\t--verbose, -v: verbose mode.")

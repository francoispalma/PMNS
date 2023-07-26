from sage.all import GF, PolynomialRing, ZZ, matrix, xgcd, gcd, factor, pari
from genamns import norm_one_of_matrix
from commonpmns import primesdict
from findm import findm

import datetime
import signal
import gc
from math import log2

Polydict = {
	"X^{n} - X^{n//2} + 1":(lambda n: 3*n/2, lambda n:[-1] + [0] * (n//2 - 1) + [1]), "X^{n} + X^{n//2} + 1":(lambda n: 3*n/2, lambda n:[-1] + [0] * (n//2 - 1) + [-1]),
	"X^{n} - X - 1":(lambda n: 2*n - 1, [1, 1]),"X^{n} - X + 1":(lambda n: 2*n - 1, [-1, 1]), "X^{n} + X - 1":(lambda n: 2*n - 1, [1, -1]), "X^{n} + X + 1":(lambda n: 2*n - 1, [-1, -1])
	}

delta = 0
phi = 64
PHI = 1<<phi
psize = 2048
primes = primesdict[psize]
deb = 0
pari.allocatemem(2**32, 2**33)

#4096128 173
#409664 43
#8192128 4
#819264 0

def sighandler(signum, frame):
	print(f"Signal called:{signum}")
	gc.collect()

signal.signal(signal.SIGTERM, sighandler)
signal.signal(signal.SIGSEGV, sighandler)

nthroot = lambda p, n: 2**(log2(p)/n)

# gamma^n = lambda
# lambda^(p-1) = 1
# lambda^(p-1)/n = gamma^(p-1) = 1
# lambda^(p-1)/n = lambda
# lambda^(p-1) = lambda^n

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
	return 4 * w * (1.02)**n * nthrootofp * ((1 + delta)**2)

def get_init_n(p, PHI, delta=0):
	init_n = p.bit_length()//int(log2(PHI)) + 1
	while PHI < get_phi_bound(init_n, init_n, nthroot(p, init_n), delta):
		init_n += 1
	if PHI < get_phi_bound(init_n, 3*init_n/2, nthroot(p, init_n), delta):
		init_n += 1
	return 42 #init_n


init_n = get_init_n(primes[0], PHI, delta)


def checkpoly(p, rts, n, w, lam, nthrootofp, delta):
	ZZX = PolynomialRing(ZZ, "X")
	E = ZZX(f"X^{n}") - ZZX(lam)
	if len(rts) == 0 or len(factor(E)) == 1:
		return []

	for rt in rts:
		gamma = int(rt[0]) % p
		if gamma <= nthrootofp or gamma >= p - nthrootofp:
			continue
		if type(lam) == int and pow(gamma, n, p) != (lam % p):
			continue

#		print("#####################################################")
#		print(datetime.datetime.now().time())
#		print("Checking for gamma =")
#		print(gamma)
#		print("#####################################################")
		B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
		print("Starting LLL:", datetime.datetime.now().time())
		Blll = list((matrix(ZZ, B).LLL()))
		print("LLL over:", datetime.datetime.now().time())
#		print(Blll[0])
#		Bbkz = list((matrix(ZZ, Blll).BKZ(algorithm="NTL", fp='rr')))
#		if type(lam) == int:
#			four_double_u = round(4*w)
#			try:
#				Mlll, M1lll = findm(p, n, gamma, lam, 64, [[int(elem) for elem in lig] for lig in Blll], B, PHI)
#			except ValueError:
#				return []
#			summaMlll = abs(Mlll[0]) + abs(lam) * sum([abs(elem) for elem in Mlll[1:]])
##			try:
##				Mbkz, M1bkz = findm(p, n, gamma, lam, 64, [[int(elem) for elem in lig] for lig in Bbkz], B, PHI)
##			except ValueError:
##				return []
##			summaMbkz = abs(Mbkz[0]) + abs(lam) * sum([abs(elem) for elem in Mbkz[1:]])
##			if summaMlll < summaMbkz:
#			summaM = summaMlll
#			M = Mlll
#			M1 = M1lll
#			else:
#				summaM = summaMbkz
#				M = Mbkz
#				M1 = M1bkz
#			if four_double_u * summaM * ((1+delta)**2) <= PHI:
#				__tmp = int(2 * summaM)
#				rho = __tmp.bit_length()
#				return [p, n, gamma, lam, rho, M, M1]
#		else:
#			return []
#		return []
#		n1Blll = int(norm_one_of_matrix(Blll))
#		n1Bbkz = int(norm_one_of_matrix(Bbkz))
#		if n1Bbkz < n1Blll:
#			B = Bbkz
#			n1B = n1Bbkz
#		else:
#			B = Blll
#			n1B = n1Blll
		B = [[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]
		B = list((matrix(ZZ, B).LLL()))
		n1B = int(norm_one_of_matrix(B))
		four_double_u = round(4*w)
		if four_double_u*n1B*((1+delta)**2) <= PHI:
			B1 = list(matrix(B).inverse() % PHI)
			__tmp = int(2 * n1B)
			rho = __tmp.bit_length()
			return [p, n, gamma, lam, rho, B, B1]
#		print("Done")
#		print(datetime.datetime.now().time())
	return []


def check_eval_strE(eval_strE, polK, p, n, w, lam, nthrootofp, delta):
	E = polK(eval_strE)
	try:
		rts = E.roots(ring=polK)
	except cypari2.handle_error.PariError:
		pari.allocatemem()
		rts = E.roots(ring=polK)
	check = checkpoly(p, rts, n, w, lam, nthrootofp, delta)
	if check:
		print(i, eval_strE)
		resdict[eval_strE] = resdict.get(eval_strE, 0) + 1
	return check

wompwomp = 0

def genpmns(p, PHI, init_n=None, amns_only=False, max_lambda=1<<64, delta=0):
	global wompwomp
	if init_n is None:
		init_n = get_init_n(p, PHI, delta)
	K = None
	polK = None
	check = []
	n = init_n
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
			check = checkpoly(p, [[int(pow(2, u, p))]], n, 2*n - 1, 2, nthrootofp, delta)
			if check:
				founddict[f"{n}"] = founddict.get(f"{n}", 0) + 1
				print(f"{i} X^{n} - 2 found pgcd 1")
				break
		u = findeasynthroot(p, n)
		if u:
			check = checkpoly(p, [[int(pow(2, u, p))]], n, 2*n - 1, 2, nthrootofp, delta)
			if check:
				founddict[f"{n}"] = founddict.get(f"{n}", 0) + 1
				print(f"{i} X^{n} - 2 found easy")
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
							check = checkpoly(p, [[int(pow(lamb, u, p))]], n, w(n), lamb, nthrootofp, delta)
							if check:
								founddict[f"{n}"] = founddict.get(f"{n}", 0) + 1
								print(f"{i} X^{n} - {lamb} found")
								break
						if check:
							break
						wompwomp += 1
				else:
					lam += 1
					continue
				break
		if check:
			break
		if K is None:
			print("Generating GF(p)", datetime.datetime.now().time())
			K = GF(p)
			print("Generating polK", datetime.datetime.now().time())
			polK = PolynomialRing(K, "X")
			print("factors", datetime.datetime.now().time())
		lam = 1
		while lam < max_lambda and PHI >= get_phi_bound(n, w(n), nthrootofp, delta):
			print("Checking lambda", lam, datetime.datetime.now().time())
			eval_strE = f"X^{n} - {lam}"
			check = check_eval_strE(eval_strE, polK, p, n, w(n), lam, nthrootofp, delta)
			if check:
				break
			eval_strE = f"X^{n} + {lam}"
			check = check_eval_strE(eval_strE, polK, p, n, w(n), -lam, nthrootofp, delta)
			if check:
				break
			lam += 1
		if check:
			break

#		print("No amns")
		#PMNS
#		for strE in polydict:
#			if "//" in strE and n & 1:
#				continue
#			w, lam = polydict[strE]
#			w = w(n)
#			if callable(lam):
#				lam = lam(n)
#			eval_strE = eval(f'f"{strE}"')
#			check = check_eval_strE(eval_strE, polK, p, n, w, lam, nthrootofp, delta)
#			if check:
#				#print(check)
#				break
		else:
			n += 1
			print("Checking n =", n, "and freeing", gc.collect())
			if n  == 43: #> 2 * init_n:
				break
			continue
		break
#	with open(f"generated/pmns{psize}{phi}{'' if delta == 0 else 'delta' + str(delta)}.py", "a") as FILE:
#		FILE.write(f"pmnsdict[{p}] = {check}\n")
	if check != []:
		print(check)
		exit()
	print("finished", datetime.datetime.now().time())
	

resdict = {}
founddict = {}

if __name__ == "__main__":
	print(psize, phi, init_n, delta, "deb =", deb)
#	if deb == 0:
#		with open(f"generated/pmns{psize}{phi}{'' if delta == 0 else 'delta' + str(delta)}.py", "w+") as FILE:
#			FILE.write("pmnsdict = {}\n")

	for i in range(deb, len(primes)):
		print()
		print(datetime.datetime.now().time(), "Checking for p =")
		print(primes[i])
		print()
		genpmns(primes[i], PHI, delta=delta)

	print(psize, phi, init_n, delta)
	print(resdict)
	print(founddict)
	print(wompwomp)

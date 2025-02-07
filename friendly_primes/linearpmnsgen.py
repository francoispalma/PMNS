#!/usr/bin/env python3
import sys

from sage.all import ZZ, next_probable_prime, is_prime, gcd, floor, ceil, log, sqrt, randrange

log2 = lambda X:log(X,2)

VERBOSE = False

def eprint(*args, **kwargs):
	if VERBOSE:
		print(*args, file=sys.stderr, **kwargs)

def gen_pmns(psize, n, alpha=1, k=1, PHI=64, max_lambda=2**64, mindelta=0, sequential=False, sparse=False, continuous=False):
	eprint(f"Finding p of size {psize} bits with n = {n} such that {str(k) + '*' if k != 1 else ''}p = {'' if alpha == 1 else str(alpha)+'*'}gamma^n - {'lambda'}")
	foundflag = False
	if PHI != 64:
		eprint(f"With phi = 2^{phi}")
	PHI = 2**PHI
	if mindelta != 0:("Accepting a minimum delta of", mindelta)
	psi = 0
	z = k
	while not(z % 2):
		psi += 1
		z//=2
	if psi:
		eprint("Warning: k being even means less potential candidates.")
		aleph = alpha//gcd(alpha,2**psi)
		psadi = psi-int(log2(gcd(alpha,2**psi)))
	else:
		aleph = alpha
		psadi = 0
	lowgamma = int(ceil((k*2**(psize-1)//alpha)**(ZZ(1)/n)))
	if sparse:
		lowgamma = ceil(lowgamma/2**32)*2**32
	assert alpha*lowgamma**n >= k*2**(psize-1)
	highgamma = int(ceil((k*2**(psize)//alpha)**(ZZ(1)/n)))
	if sparse:
		highgamma = ceil(highgamma/2**32)*2**32
	assert alpha*highgamma**n >= k*2**(psize)
	eprint(f"{(highgamma - lowgamma)//(2**32 if sparse else 1)} potential values for gamma")
	if lowgamma == highgamma:
		return 0
	if alpha == 1:
		lambda_max = lambda gamma : (floor(sqrt((((n-1)*(gamma-2))**2 - 2*(n-1)*(gamma - PHI//((1+mindelta)**2) - 1) + 1)) - ((n-1)*(gamma-2)+1)))//(2*(n-1))
	else:
		lambda_max = lambda gamma : max(alpha,floor( (PHI/(2*(aleph*gamma - 1)*(1+mindelta)**2) - alpha)/(n-1)))
	if lambda_max(lowgamma) < 1:
		eprint("lambda_max for low gamma < 1")
		return 0
	eprint(f"Maximum potential lambda of {lambda_max(lowgamma)}")
	try:
		avglambda = round(k*log(1.5*2**(psize-1))/2)
		eprint(f"Average expected value for lambda of {avglambda}")
		expected_value = lambda_max((highgamma+lowgamma)//2)/avglambda/2**psadi
		expected_str = float(expected_value)
		if expected_str > 1:
			expected_str = round(expected_str, 4)
		eprint(f"Each gamma has an expected number of valid primes of {expected_str}")
		expected_str = float(expected_value*(highgamma-lowgamma))
		if expected_str > 1:
			expected_str = round(expected_str, 4)
		eprint(f"Expected total number of primes in the interval is {expected_str}")
	except OverflowError:
		pass
	gamma = lowgamma - (2**32 if sparse else 1)
	display = 2**64
	eprint("\nOutput format is (p, n, gamma, rho, E, delta)")
	while True:
		eprint(f"\r{gamma - lowgamma}\t{display}                       ", end="\r")
		if sequential:
			gamma += 2**32 if sparse else 1
			if gamma == highgamma:
				break
		else:
			gamma = randrange(lowgamma, highgamma, 2**32 if sparse else 1)
		lambdamax = min(lambda_max(gamma), max_lambda)
		if lambdamax >= 1:
			# We use probable for speed, we check with is_prime later.
			p = next_probable_prime(floor(ZZ(alpha*gamma**n - lambdamax)/k))
			lambda_ = alpha*gamma**n - k*p
			# We try to find the smallest valid lambda.
			prev_l = lambda_
			prev_p = p
			pots = [lambda_]
			while lambda_ > 0:
				prev_p = p
				p = next_probable_prime(p)
				prev_l = lambda_
				lambda_ = alpha*gamma**n - k*p
				pots += [lambda_]
			if abs(prev_l) < abs(lambda_):
				lambda_ = prev_l
				p = prev_p
			if abs(lambda_) < display:
				display = abs(lambda_)
			if not continuous:
				pots = [lambda_]
			for lambda_ in pots:
				lamed = abs(lambda_)//(2**psi)
				w = max(alpha*n, alpha + abs(lambda_)*(n-1))
				normoneofG = max((alpha*gamma)//(2**psi) + 1, gamma + lamed)
				delta = floor(sqrt(ZZ(PHI)/(2*w*(normoneofG-2))))-1
				if abs(lambda_) <= lambdamax and gcd(lambda_, 2**psi) == 2**psi and gcd(alpha*gamma, 2**psi) == 2**psi and delta >= mindelta and is_prime(p):
					foundflag = True
					rho = normoneofG - 1
					E = f"{'' if alpha == 1 else str(alpha)+'*'}X^{n} {'-' if lambda_>0 else '+'} {abs(lambda_)}"
					if not sparse or gcd(alpha*gamma//2**psi, 2**32) == 2**32:
						if not continuous:
							return (p, n, gamma, rho, E, delta)
						else:
							print((p, n, gamma, rho, E, delta))
	if foundflag:
		return "Over."

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
	opts = [arg for arg in args if arg.startswith("-")]
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
			assert max_lambda >= 1, f"Invalid value for max_lambda: {max_lambda}"
			max_lambda = min(2**phi, max_lambda)
		except Exception as e:
			print(e)
			exit()
	else:
		max_lambda = 2**64
	SEQUENTIAL = False
	if "-s" in args:
		SEQUENTIAL = True
		args.pop(args.index("-s"))
	if "--sequential" in args:
		SEQUENTIAL = True
		args.pop(args.index("--sequential"))
	SPARSE = False
	if "-S" in args:
		SPARSE = True
		args.pop(args.index("-S"))
	if "--sparse" in args:
		SPARSE = True
		args.pop(args.index("--sparse"))
	CONTINUOUS = False
	if "-c" in args:
		CONTINUOUS = True
		args.pop(args.index("-c"))
	if "--continuous" in args:
		CONTINUOUS = True
		args.pop(args.index("--continuous"))
	arguments = [arg for arg in args if not arg.startswith("-")]
	if len(arguments) >= 2:
		try:
			psize = int(arguments[1])
		except ValueError:
			print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
			exit()
		if psize <= phi:
			print("Prime Size too small")
			exit()
		if len(arguments) >= 3:
			try:
				n = int(arguments[2])
			except ValueError:
				print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
				exit()
			if n <= 1:
				print("Invalid value for n:", n)
				exit()
		else:
			n = psize//phi
			if SPARSE:
				w = lambda n:3*n - 2
			else:
				w = lambda n:2*n - 1
			while 2**phi <= 2*w(n)*(ceil(2**((psize-1)/ZZ(n))))*((1+delta)**2):
				n += 1
				if n > psize:
					break
		alpha = 1
		if len(arguments) >= 4:
			try:
				alpha = int(arguments[3])
			except ValueError:
				print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
				exit()
			if alpha == 0:
				print("Invalid value for alpha: 0")
				exit()
		k = 1
		if len(arguments) >= 5:
			try:
				k = int(arguments[4])
			except ValueError:
				print("Invalid arguments:", arguments[1:].__repr__().replace("'", "")[1:-1])
				exit()
			if k == 0:
				print("Invalid value for k: 0")
				exit()
		check = gen_pmns(psize, n, alpha, k, phi, max_lambda, delta, SEQUENTIAL, SPARSE, CONTINUOUS)
		if check:
			print(check)
		else:
			print("No results found with these parameters")
	else:
		print("Usage: ./lineargenpmns.py ",end="")
		underline = lambda X: '\033[4m' + X + '\033[0m'
		print(underline("PRIMESIZE"), end=" [")
		print(underline("N"), end="] [")
		print(underline("ALPHA"), end="] [")
		print(underline("K"), end="] [")
		print(underline("OPTION"), end="]...\n\n")
		print("\t--delta, -d: use specified delta for number of free additions.")
		print("\t--phi, -F: sets the parameter PHI = 2^phi. Default value is phi=64.")
		print("\t--lambda, -l: use specified value for maximum lambda allowed.")
		print("\t--sequential, -s: iterate over all gamma in order and not randomly.")
		print("\t--sparse, -S: only allow DoubleSparse outputs.")
		print("\t--continuous, -c: keep generation going forever (or until finished).")
		print("\t--verbose, -v: verbose mode.")

import sys
import os

from ast import literal_eval
from contextlib import redirect_stdout

from commonpmns import pmnsWBdicts, primesdict
from precalcs import do_precalcs
from precalcs128 import do_precalcs as do_precalcs128
from ops import horner_modulo

def do_checks(Psize="1024", phi=""):
	# Parameters
	nbmults = 2**20
	pmnsdict = pmnsWBdicts[int(Psize + phi)]
	primes = primesdict[int(Psize)]
	precalc = do_precalcs128 if phi == "128" else do_precalcs
	filename = "results/checks" + Psize + phi
	sphi = "p128.exe" if phi == "128" else "main.exe"
	PHI = 2**128 if phi == "128" else 2**64

	# We wipe the files
	with open(filename, "w+") as f:
		pass

	# We run the checks
	for i in range(len(pmnsdict)):
		(p, n, gamma, lam, rho, M, M1) = pmnsdict[primes[i]]
		RHO = 2**rho
		with open("params" + phi + ".h", "w+") as f:
			with redirect_stdout(f):
				precalc(p, n, gamma, lam, rho, M, M1)
		os.system(f"make {sphi}")
		os.system(f"./{sphi} {nbmults} > log{phi}")
		counter = 0
		mtr = 0
		with open(f"log{phi}", "r") as f:
			while True:
				try:
					a = literal_eval(f.readline()[:-1])
					b = literal_eval(f.readline()[:-1])
					c = literal_eval(f.readline()[:-1])
					c_check = (horner_modulo(a, gamma, p) * horner_modulo(b, gamma, p) * pow(PHI, -1, p)) % p
					if horner_modulo(c, gamma, p) != c_check:
						counter += 1
					for elem in c:
						if abs(elem) >= RHO:
							mtr += 1
#							print("More than Rho:", abs(elem), RHO)
#							exit()
				except SyntaxError:
					break
		os.system(f"echo \"({counter}, {mtr})\" >> {filename}")

if __name__ == "__main__":
	if len(sys.argv) >= 2:
		try:
			psize = int(sys.argv[1])
			if psize not in list(primesdict.keys()):
				print("Prime Size not handled")
				exit()
			try:
				phi = sys.argv[2]
			except IndexError:
				phi = ""
			if phi not in ["128", "64", ""]:
				print("Value of Phi not handled")
				exit()
			elif phi == "64":
				phi = ""
		except ValueError:
			print("Invalid arguments")
			exit()
		do_checks(str(psize), phi)
	else:
		print("Not enough arguments: Psize [Phi]")

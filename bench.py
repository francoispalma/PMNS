import sys
import os

from contextlib import redirect_stdout

from commonpmns import pmnsdicts, primesdict
from precalcs import do_precalcs
from precalcs128 import do_precalcs as do_precalcs128

def do_bench(Psize="1024", phi=""):
	# Parameters
	pmnsdict = pmnsdicts[int(Psize + phi)]
	primes = primesdict[int(Psize)]
	precalc = do_precalcs128 if phi == "128" else do_precalcs
	filename = "results/results" + Psize + phi

	# We wipe the files
	with open(filename, "w+") as f:
		pass
	with open(filename + "pre", "w+") as f:
		pass
	if phi == "128":
		with open(filename + "hyb", "w+") as f:
			pass

	# We run the benchmark
	for i in range(len(pmnsdict)):
		(p, n, gamma, lam, rho, M, M1) = pmnsdict[primes[i]]
		with open("params" + phi + ".h", "w+") as f:
			with redirect_stdout(f):
				precalc(p, n, gamma, lam, rho, M, M1)
		os.system("make bench" + phi + ".exe")
		os.system("./bench" + phi + ".exe >> " + filename)
		os.system("./bench" + phi + ".exe pre >> " + filename + "pre")
		if phi == "128":
			os.system("./bench" + phi + ".exe hyb >> " + filename + "hyb")

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
			do_bench(str(psize), phi)
		except ValueError:
			print("Invalid arguments")
	else:
		print("Not enough arguments: Psize [Phi]")

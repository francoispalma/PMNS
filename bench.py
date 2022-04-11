import os

from contextlib import redirect_stdout

from generatedpmns import pmnsdict
from primes1024 import PRIMES1024 as primes
from precalcs import do_precalcs

with open("results", "w+") as f:
	pass
for i in range(100):
	(p, n, gamma, lam) = pmnsdict[primes[i]]
	with open("params.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs(p, n, gamma, lam)
	os.system("make")
	os.system("./main.exe >> results")

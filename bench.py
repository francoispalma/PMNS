import os

from contextlib import redirect_stdout

from generated4096pmns import pmnsdict
from primes4096 import PRIMES4096 as primes
from precalcs import do_precalcs

with open("results4096", "w+") as f:
	pass
for i in range(len(pmnsdict)):
	(p, n, gamma, lam) = pmnsdict[primes[i]]
	with open("params.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs(p, n, gamma, lam)
	os.system("make bench.exe")
	os.system("./bench.exe >> results4096")

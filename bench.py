import os

from contextlib import redirect_stdout

from generated2048pmns import pmnsdict
from primes2048 import PRIMES2048 as primes
from precalcs import do_precalcs

with open("results2048", "w+") as f:
	pass
for i in range(1000):
	(p, n, gamma, lam) = pmnsdict[primes[i]]
	with open("params.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs(p, n, gamma, lam)
	os.system("make bench.exe")
	os.system("./bench.exe >> results2048")

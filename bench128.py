import os

from contextlib import redirect_stdout

from generated2048pmns128 import pmns128dict
from primes2048 import PRIMES2048 as primes
from precalcs128 import do_precalcs as do_precalcs128

with open("results2048128", "w+") as f:
	pass
for i in range(1000):
	(p, n, gamma, lam) = pmns128dict[primes[i]]
	with open("params128.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs128(p, n, gamma, lam)
	os.system("make bench128.exe")
	os.system("./bench128.exe >> results2048128")

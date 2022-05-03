import os

from contextlib import redirect_stdout

from generated4096pmns128 import pmns128dict
from primes4096 import PRIMES4096 as primes
from precalcs128 import do_precalcs as do_precalcs128

with open("results4096128", "w+") as f:
	pass
for i in range(10):
	(p, n, gamma, lam) = pmns128dict[primes[i]]
	with open("params128.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs128(p, n, gamma, lam)
	os.system("make bench128.exe")
	os.system("./bench128.exe >> results4096128")

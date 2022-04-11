import os

from contextlib import redirect_stdout
from sage.all import matrix, ZZ
from sage.modules.free_module_integer import IntegerLattice
from math import ceil

#from generatedpmns import pmnsdict, pmns128dict
from generatedpmns128 import pmns128dict
from primes1024 import PRIMES1024 as primes
from precalcs128 import do_precalcs as do_precalcs128
from precalcs import do_precalcs

with open("results128", "w+") as f:
	pass
for i in range(1000):
	(p, n, gamma, lam) = pmns128dict[primes[i]]
	with open("params128.h", "w+") as f:
		with redirect_stdout(f):
			do_precalcs128(p, n, gamma, lam)
	os.system("make p128.exe")
	os.system("./p128.exe >> results128")

import sys, os

from sage.all import is_prime
from math import log2, ceil
from time import perf_counter
from contextlib import redirect_stdout

from genamns import gen_amns
from precalcs128 import do_precalcs as precalc128
from precalcs import do_precalcs as precalc


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("Expected format: python3 completegen.py PRIME [PHI]")
		print("With PRIME the prime integer to be used for the modular reductions.")
		print("With PHI the size of registers you want used (64 or 128, default is 64).")
	else:
		try:
			p = int(sys.argv[1])
		except ValueError:
			print(f"Error: expected a prime number but got '{sys.argv[1]}' as argument instead.")
			exit()
		if not is_prime(p):
			print(f"Error: integer {p} is not prime.")
			exit()
		phi = 64
		if len(sys.argv) > 2:
			try:
				phi = int(sys.argv[2])
			except ValueError:
				phi = 64
		if phi not in [64, 128]:
			phi = 64
		print(f"Generating PMNS for integer of size {ceil(log2(p))}...\n")
		STARTSTAMP = perf_counter()
		pmns = gen_amns(p, phi)
		STOPSTAMP = perf_counter()
		print(f"Generation done in {str(STOPSTAMP - STARTSTAMP)[:6]} seconds.")
		p, n, gamma, lam, rho, M, M1 = pmns
		print(f"The generated pmns has polynomials of degree {n}")
		if not os.path.exists("output"):
			os.makedirs("output")
		sphi = "128" if phi == 128 else ""
		with open(f"output/params{sphi}.h", "w+") as f:
			with redirect_stdout(f):
				if phi == 64:
					precalc(p, n, gamma, lam, rho, M, M1)
				elif phi == 128:
					precalc128(p, n, gamma, lam, rho, M, M1)
		os.system(f"cp pmns{sphi}.[ch] output/")
		os.system("cp utilitymp.[ch] output/")
		os.system("cp structs.[ch] output/")
		with open("output/makefile", "w+") as f:
			with redirect_stdout(f):
				print(f"""FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -Wno-restrict -lgmp
CC = gcc

all: main.exe

main.exe: main.c pmns{sphi}.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns{sphi}.o: pmns{sphi}.c pmns{sphi}.h
	$(CC) -c $< $(FLAGS)

structs.o: structs.c structs.h
	$(CC) -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

demo: main.exe
	./main.exe 41781314183145 734865131442

""")
		with open("output/main.c", "w+") as f:
			with redirect_stdout(f):
				print("""#include <stdio.h>

int main(int argc, char** argv)
{
	if(argc > 2)
		printf("%s %s\\n", argv[1], argv[2]);
	else
		printf("To output a^b %% p please call this executable with the numbers a and b as parameters\\n");
	return 0;
}
""")
		print("\nGenerated code can be found in the directory output.")
		print("A demonstration can be launched with the command 'make demo'")
		print("The demo executable outputs arg1^arg2 % p")
		
# 101446619888572739133666196758813369433960395136978547326876308901415794561163

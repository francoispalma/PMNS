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
		print(f"Generating PMNS {phi} bits for integer of size {ceil(log2(p))}...\n")
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
				print(f"""FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -Wno-restrict -Wno-unused-variable -lgmp
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
	./main.exe 0x417a131f18c14b 0xa348e5131d42

""")
		with open("output/main.c", "w+") as f:
			with redirect_stdout(f):
				print(f"""#include <stdio.h>
#include \"pmns{sphi}.h\"
#include \"utilitymp.h\"

int main(int argc, char** argv)
{{
	if(argc > 2)
	{{
		poly{sphi} A, C;
		poly B, aux;
		init_polys(1, &B, &aux, NULL);
		init_poly{sphi}s(N, &A, &C, NULL);
		
		convert_string_to_poly(&B, argv[2]);
		convert_string_to_amns{sphi}(A, argv[1]);
		
		amns{sphi}_montg_ladder(C, A, B);
		
		convert_amns{sphi}_to_poly(&aux, C);
		
		printf("0x");
		mp_print(aux);
		
		free_polys(B, aux, NULL);
		free_poly{sphi}s(A, C, NULL);""")
				print("""\t}
	else
		printf("To output a^b %% p please call this executable with the numbers a and b as parameters\\n");
	return 0;
}
""")
		print("\nGenerated code can be found in the directory output.")
		print("A demonstration can be launched with the command 'make demo'")
		print("The demo executable outputs arg1^arg2 % p (arguments in hex)")
		
# 101446619888572739133666196758813369433960395136978547326876308901415794561163

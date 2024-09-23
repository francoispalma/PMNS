# PMNS using friendly primes
This repository contains codes to generate Polynomial Modular Number Systems (PMNS) for a friendly class of primes and code to generate parameters for use in the accompanying C code for it. It also contains benchmarking code to compare the performances with other modular arithmetic of a given size in bits.

The subdirectory 'curves' contains example of Elliptic curves generated from some randomly generated primes belonging to the PMNS-friendly class.

To generate a table from the paper simply enter:
> make tableX

where X is the table number. Note that you need to disable turbo boost to get accurate results. To disable turbo boost on linux with an intel processor, one can do the following:
> sudo ./disableturbo.sh

The PMNS generator will display the expected syntax when executed like so:
> ./linearpmnsgen.py

The code generator expects either a PMNS as parameter or a filename and line number if using a file in which many PMNS are stored. The format expected is the one from the output of linearpmnsgen, that is to say:
> ./codegen.py p, n, gamma, rho, E, delta

For example:
> ./codegen.py 57896044618658097711785492504343953926634992332820282019728792003956564819949, 5, 2251799813685248, 2251799813685266, 'X^5 - 19', 0

or alternatively:
> ./codegen.py [filename] [linenumber]

will generate the appropriate code into the file "params.h" which is used by the "hpmns.c" file, which is intended to be imported in other projects as opposed to "pmns.c" which has been specifically written for the benchmarks of the paper. An example of main file using hpmns.c can be seen in the file "hmain.c" and an example compilation can be seen by performing:
> make hmain.exe

which can then be executed with the line
> ./hmain.exe

to benchmark the resulting PMNS and check the result of a million multiplications.

Note : To use the PMNS generator, you will need SageMath library which can be found here: http://www.sagemath.org/

# Curves constructed on PMNS-friendly primes
This subdirectory contains example of Elliptic curves generated from some randomly generated primes belonging to the PMNS-friendly class.

Each file of the name curves[size].py contains example of curves constructed on prime fields of the size given in the file name. The curve parameters are arranged in python dictionaries for ease of use.

The format used is curvedict[prime] = [b, trace, embedding degree].

All curves are of prime order using the short weierstrass equation y² = x³ + ax + b with a = -3 and verify the brainpool criteria for Elliptic Curve Cryptography.

Each file of the name pmns[size].py contains the accompanying PMNS constructed on the primes used for the curve generation. For each prime used as a key in a curvedict in curves[size].py, there is a corresponding pmns in pmns[size].py arranged in a python dictionary for ease of use.

The format used is pmnsdict[prime] = (p, n, gamma, rho, E, delta).

The format is the same as the output from linearpmnsgen.py from the parent directory, and the same as the one expected as parameter by codegen.py in the parent directory. For example one can run
> ./codegen.py curves/pmns256.py 2

from the parent directory to generate code for the first PMNS corresponding to the first curve in curves256.py.

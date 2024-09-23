# Curves constructed on PMNS-friendly primes
This subdirectory contains example of Elliptic curves generated from some randomly generated primes belonging to the PMNS-friendly class.

Each file of the name curves[size].py contains example of curves constructed on prime fields of the size given in the file name. The curve parameters are arranged in python dictionaries for ease of use.

The format used is curvedict[prime] = [b, trace, embedded degree].

All curves are of prime order using the short weierstrass equation y² = x³ + ax + b with a = -3 and verify the brainpool criteria for Elliptic Curve Cryptography.

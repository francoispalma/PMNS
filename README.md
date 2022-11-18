# PMNS

This project is an implementation of the Polynomial Modular Number System (PMNS), a system for fast modular arithmetic operations for use in cryptography.

Focus is put on the Adapted Modular Number System (AMNS) so far. Generation is done for those for now. Support for other PMNS exists currently but accompanying generation doesn't exist as of yet.

The implementation itself is in C with utility functions in Python for generation purposes.

Currently a 64 bit and 128 bit versions are available. A 104 bit version using AVX512IFMA is possible in the future.

## Requirements

Installing sagemath is required to use the code generation.
> sudo apt install sagemath-common

then

> sudo apt install sagemath

## Usage

To generate all the code with 64-bit or 128-bit implementation needed for operations with a specific prime p along with a demo code to calculate a^b % p and accompanying makefile:
> python3 completegen.py {p} [PHI=64 or 128]

To load a specific PMNS of size PSIZE and index INDEX:
> make loadpmns [PSIZE=size] [INDEX=index]

To load a specific PMNS of size PSIZE and index INDEX using basis reduction:
> make loadpmnsWB [PSIZE=size] [INDEX=index]

To display a benchmark of the current version with the currently loaded PMNS:
> make bench

Alternatively for the 128 bit version:
> make bench128

For details as to the generated PMNS degrees and lambda values:
> python3 deghisto.py

For a table of current results summarized:
> python3 showresults.py


### Additional utilities for coding purposes

To check if the current version leaks memory:
> make check

To test out if the current code properly computes A times B mod P (counts the number of errors):
> make proof

Alternatively for the 128 bit version:
> make proof128

## Acknowledgements
Todo

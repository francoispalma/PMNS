# PMNS

This project is an implementation of the Polynomial Modular Number System (PMNS), a system for fast modular arithmetic operations for use in cryptography.

Focus is put on the Adapted Modular Number System (AMNS) so far but this may change in the future.

The implementation itself is in C with utility functions in Python for generation purposes.

Currently a 64 bit and 128 bit versions are available. A 104 bit version using AVX512IFMA is possible in the future.

## Usage

To display various graphics and results in a visual format:
> python3 histo.py

For details as to the generated PMNS degrees and lambda values:
> python3 deghisto.py

To load a specific PMNS of size PSIZE and index INDEX:
> make loadpmns [PSIZE=size] [INDEX=index]

To display a benchmark of the current version with the currently loaded PMNS:
> make bench

Alternatively for the 128 bit version:
> make bench128

### Additional utilities for coding purposes

To check if the current version leaks memory:
> make check

To test out if the current code properly computes A times B mod P (counts the number of errors):
> make proof

Alternatively for the 128 bit version:
> make proof128

## Acknowledgements
Todo

# Multi-precision PMNS with CIOS reduction

This folder contains the code used to measure cycle counts in the relevant tables in the "Multi-precision PMNS with CIOS reduction" paper. You will need GCC and GnuMP to compile the codes. For the performances to be accurate on intel processors the Turbo Boost® needs to be disabled. For linux we provide the script "disableturboboost.sh" to disable the turbo boost. To run it:
> chmod +x disableturboboost.sh && sudo ./disableturboboost.sh

## Benchmarks
To run the benchmark for table 5:
> make table5

Similarly, to run the benchmark for table 6 (your processor needs the cpuflag avx512ifma or the code won't run):
> make table6

Similarly, to run the benchmark for table 8 (your processor needs the cpuflag avx512ifma or the code won't run):
> make table8

Similarly, to run the benchmark for table 9 (your processor needs the cpuflag avx512ifma or the code won't run):
> make table9

## Generation
We also make a multi-precision PMNS generation code available as well as an accompanying parameter generation script. To use them, you will need the SageMath library which can be found here: http://www.sagemath.org/

To use the generation code simply edit the values at the top of the python script (genmppmns.py) as needed. Then simply paste the output PMNS (p,n,gamma,lambda,rho) in (avx512)codegen.py at the appropriate spot. Run the generation script to generate a .h parameter file for the appropriate .c file.


## Correctness and bounds check
This repository also contains code to check that each PMNS is valid and has correct parameter consistency. The code relies on the SageMath library which can be found here: http://www.sagemath.org/

The checks are split by table in the makefile but the code is a simple python script that can be reused as needed.

To check that the PMNS used in table X are valid and that operations don't produce results whose coefficients are more than rho simply run:
> make checktableX

For example for table 5:
> make checktable5

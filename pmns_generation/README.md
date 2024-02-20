# PMNS generation

This folder contains the code to generate PMNS.

## Requirements

Installing sagemath is required to use the code generation.
> sudo apt install sagemath-common

then

> sudo apt install sagemath

GCC is needed to compile the codes (although one may edit the resulting makefiles for any other compiler, other compilers have not been tested and may not function properly).

## Instructions

To generate a PMNS for a specific prime:

> ./gen_pmns.py -p {prime}

To generate a PMNS for a specific prime size in bits:

> ./gen_pmns.py {prime size}

More options available such as verbose mode and various generation parameters.

For a list of options:

> ./gen_pmns.py

Once a PMNS is generated, generating the code is done as such:

### Generating pmns from a file

You can redirect the output from the generation into a file and then generate with that file as an input as follows:

> ./gen_pmns.py [OPTIONS] > tmp
> ./codegen.py tmp

After the codegen is run a params.h file is created with all the necessary parameters to compile pmns.c

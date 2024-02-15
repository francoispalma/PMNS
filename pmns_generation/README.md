# PMNS generation

This folder contains the code to generate PMNS. You will need Sagemath to run the generation (see parent folder for instructions on how to install Sagemath) and GCC to compile the codes (although one may edit the resulting makefiles for any other compiler, other compilers have not been tested and may not function properly).

## Instructions

To generate a PMNS for a specific prime:

> ./gen_pmns.py -p {prime}

To generate a PMNS for a specific prime size in bits:

> ./gen_pmns.py {prime size}

More options available such as verbose mode and various generation parameters.

For a list of options:

> ./gen_pmns.py

Once a PMNS is generated, generating the code is done as such:

> TODO

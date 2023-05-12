# PMNS for cryptography: a guided tour

This folder contains the code to generate each of the tables in the "PMNS for cryptography: a guided tour" paper. You will need GCC and GnuMP to compile the codes.

For each folder to compile and run the code to generate the table:
make && ./main.exe

Table 4 folder includes code to check the modular multiplication in a PMNS to show that we are computing it correctly.

Table 6 folder includes code to show the parallel toeplitz recursion method is also correct in its computations.

Table 9 folder includes code to compare the result after a modular multiplication using the toeplitz recursive method in all three operations of the PMNS vs the GnuMP result and thus shows the result is correct.

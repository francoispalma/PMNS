CC = gcc-12

redintcomp.exe: redintcomp.c redparams.h
	$(CC) -o redintcomp.exe redintcomp.c -g -Wall -Wextra -Wno-comment -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -lgmp -fopenmp

redparams.h:
	python3 codegen.py

check:
	python3 redcheck.py

CC = gcc-12
FASTFLAGS = -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize
WARNINGS = -g -Wall -Wextra -Wno-comment -Wno-restrict

equalitytest.exe: equalitytest.c
	@$(CC) -o $@ $< $(WARNINGS) $(FASTFLAGS) -fwhole-program -lgmp

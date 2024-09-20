FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize
CC = gcc

Npmns.o: pmns.c
	$(CC) -c $< $(FLAGS) -D VARIOUSNVALUES -o $@

Lpmns.o: pmns.c
	$(CC) -c $< $(FLAGS) -D LARGEPMNS -o $@

table3.exe: table3.c structs.o Npmns.o bench.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

table6.exe: table6.c structs.o Lpmns.o bench.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

table%.exe: table%.c structs.o pmns.o eccoptimizedcode.o bench.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

hpmns.o: hpmns.c params.h
	$(CC) -c $< $(FLAGS)

%.o: %.c
	$(CC) -c $< $(FLAGS) -Wno-unused-function

table%: table%.exe
	./$<

hmain.exe: hmain.c structs.o hpmns.o bench.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

conversioncheck.exe: conversioncheck.c structs.o hpmns.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

clean:
	rm -rf *.o
	rm -rf *.exe

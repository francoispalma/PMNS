FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program
CC = gcc-12

main:
	$(CC) -o mppmns.exe mppmns.c $(FLAGS) && ./mppmns.exe
	python3 mpcheck.py

FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops -Wno-shift-negative-value -fopenmp
CC = gcc

main.exe: intel-measurement.o pmns8192.o structs.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

%.o: %.c
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

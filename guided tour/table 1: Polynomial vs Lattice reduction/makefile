FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops
CC = gcc

main.exe: intel-measurement.o pmns256.o pmns512.o pmns1024.o structs.o
	$(CC) -o $@ $^ $(FLAGS)

%.o: %.c
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

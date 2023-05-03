FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops -Wno-shift-negative-value
CC = gcc

main.exe: intel-measurement.o pmns256.o pmns512.o pmns1024.o pmns2048.o pmns4096.o pmns8192.o pmns256128.o pmns512128.o pmns1024128.o pmns2048128.o pmns4096128.o pmns8192128.o structs.o
	$(CC) -o $@ $^ $(FLAGS)

%.o: %.c
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

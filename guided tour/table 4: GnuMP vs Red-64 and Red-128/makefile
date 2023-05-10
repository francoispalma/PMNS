FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops -Wno-shift-negative-value
CC = gcc

main.exe: intel-measurement.o pmns1024.o pmns2048.o pmns4096.o pmns8192.o pmns8192128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

%.o: %.c
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

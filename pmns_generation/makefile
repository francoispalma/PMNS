FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -funroll-loops
CC = gcc

demo.exe: demo.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns.o: pmns.c pmns.h params.h
	$(CC) -c $< $(FLAGS) 

structs.o: structs.c structs.h
	$(CC) -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	$(CC) -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

demo: demo.exe
	./demo.exe

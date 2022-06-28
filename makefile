FLAGS= -Wall -Wextra -g -O3 -funswitch-loops
CC = gcc
PSIZE = 1024
INDEX = 1

all: main.exe

main.exe: main.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns.o: pmns.c pmns.h params.h
	$(CC) -c $< $(FLAGS)

params.h: precalcs.py pyparams.py
	python3 $< > $@ $(PSIZE) $(INDEX)

structs.o: structs.c structs.h
	$(CC) -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	$(CC) -c $< $(FLAGS)

hardcode.exe: hardcode.c
	$(CC) -o $@ $^ $(FLAGS)

p128.exe: main128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

pmns128.o: pmns128.c pmns128.h params128.h
	$(CC) -c $< $(FLAGS) -lgmp

params128.h: precalcs128.py pyparams128.py
	python3 $< > $@ $(PSIZE) $(INDEX)

bench.exe: intel-measurement.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

bench128.exe: intel-measurement128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

multbench.exe: multmeasurement.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

multbench128.exe: multmeasurement128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

clean:
	rm -rf *.o
	rm -rf *.exe

gedit: clean
	gedit makefile &
	gedit *.* &

check: main.exe
	valgrind -s --leak-check=full --track-origins=yes ./main.exe

demo: main.exe
	./main.exe

proof: main.exe
	./main.exe > log
	python3 proof.py

proof128: p128.exe
	./p128.exe > log
	python3 check128.py

p128: p128.exe
	./p128.exe

bench: bench.exe
	./bench.exe

prebench: bench.exe
	./bench.exe pre

bench128: bench128.exe
	./bench128.exe

prebench128: bench128.exe
	./bench128.exe pre

progress: *.c *.py *.h makefile
	git add .
	git commit -m "progress"
	git push

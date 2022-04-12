FLAGS= -Wall -Wextra -g -O3

all: main.exe

main.exe: main.c pmns.o structs.o utilitymp.o
	gcc -o $@ $^ $(FLAGS)

pmns.o: pmns.c pmns.h params.h
	gcc -c $< $(FLAGS)

params.h: precalcs.py
	python3 $< > $@

structs.o: structs.c structs.h
	gcc -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	gcc -c $< $(FLAGS)

hardcode.exe: hardcode.c
	gcc -o $@ $^ $(FLAGS)

p128.exe: main128.c pmns128.o structs.o utilitymp.o
	gcc -o $@ $^ $(FLAGS)

pmns128.o: pmns128.c pmns128.h params128.h
	gcc -c $< $(FLAGS)

params128.h: precalcs128.py pyparams.py
	python3 $< > $@

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

progress: *.c *.py *.h makefile
	git add .
	git commit -m "progress"
	git push

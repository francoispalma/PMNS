FLAGS= -Wall -Wextra -g -O3

all: main.exe

main.exe: montgom.o
	gcc -o $@ $^ $(FLAGS)

montgom.o: montgom.c montgom.h
	gcc -c $< $(FLAGS)

hardcode.exe: hardcode.o	
	gcc -o $@ $^ $(FLAGS)

hardcode.o: hardcode/hardcode.c
	gcc -c $< $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe

gedit: clean
	gedit makefile &
	gedit *.* &

check: main.exe
	valgrind -s --leak-check=full ./main.exe

demo: main.exe
	./main.exe

proof: main.exe
	./main.exe > log
	python3 proof.py

hard: hardcode.exe
	./hardcode.exe

hardproof: hardcode.exe
	./hardcode.exe > log
	python3 proof.py

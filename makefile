FLAGS= -Wall -Wextra -g -O3

all: main.exe

main.exe: main.c montgom.o structs.o utilitymp.o utilitymp_core.o
	gcc -o $@ $^ $(FLAGS)

structs.o: structs.c structs.h
	gcc -c $< $(FLAGS)

montgom.o: montgom.c montgom.h params.h
	gcc -c $< $(FLAGS)

params.h: precalcs.py
	python3 $< > $@

hardcode.exe: hardcode.o
	gcc -o $@ $^ $(FLAGS)

hardcode.o: hardcode/hardcode.c
	gcc -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	gcc -c $< $(FLAGS)

utilitymp_core.o: utilitymp_core.c utilitymp_core.h
	gcc -c $< $(FLAGS)

mppmns.o: mppmns.c mppmns.h
	gcc -c $< $(FLAGS)

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

hard: hardcode.exe
	./hardcode.exe

hardproof: hardcode.exe
	./hardcode.exe > log
	python3 proof.py

params128.h: genamns128.py
	python3 genamns128.py > params128.h

p128.exe: mppmns.o structs.o utilitymp_core.o
	gcc -o $@ $^ $(FLAGS)

p128: p128.exe
	./p128.exe

progress: *.c *.py *.h makefile
	git add .
	git commit -m "progress"
	git push

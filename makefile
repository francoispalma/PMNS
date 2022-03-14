FLAGS= -Wall -Wextra -g -O3

all: main.exe

main.exe: main.c montgom.o mppmns.o utilitymp.o
	gcc -o $@ $^ $(FLAGS)

montgom.o: montgom.c montgom.h params.h
	gcc -c $< $(FLAGS)

hardcode.exe: hardcode.o	
	gcc -o $@ $^ $(FLAGS)

hardcode.o: hardcode/hardcode.c
	gcc -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
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

tmp.txt: precalcs.py
	python3 precalcs.py > tmp.txt

tmp: tmp.txt
	make -B
	./main.exe

proof: main.exe
	./main.exe > log
	python3 proof.py

hard: hardcode.exe
	./hardcode.exe

hardproof: hardcode.exe
	./hardcode.exe > log
	python3 proof.py

mult: mult.exe
	./mult.exe

progress: *.c *.py *.h makefile
	git add .
	git commit -m "progress"
	git push

FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -Wno-restrict -funroll-loops
CC = gcc
PSIZE = 1024
INDEX = 0

all: main.exe

main.exe: main.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns.o: pmns.c pmns.h params.h
	$(CC) -c $< $(FLAGS) 

params.h: precalcs.py
	python3 $< > $@

structs.o: structs.c structs.h
	$(CC) -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	$(CC) -c $< $(FLAGS)

p128.exe: main128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns128.o: pmns128.c pmns128.h params128.h
	$(CC) -c $< $(FLAGS)

params128.h: precalcs128.py
	python3 $< > $@

bench.exe: intel-measurement.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

bench128.exe: intel-measurement128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

multbench.exe: multmeasurement.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

multbench128.exe: multmeasurement128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

expbench.exe: modexpmeasures.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

expbench128.exe: modexpmeasures128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe
	rm -rf output/

gedit: clean
	gedit makefile &
	gedit *.* &

check: main.exe
	valgrind -s --leak-check=full --track-origins=yes ./main.exe 100

demo: main.exe
	./main.exe 1

proof: main.exe
	./main.exe 100 > log
	python3 proof.py

proof128: p128.exe
	./p128.exe 100 > log128
	python3 check128.py

p128: p128.exe
	./p128.exe 100

bench: bench.exe
	./bench.exe

bench128: bench128.exe
	./bench128.exe

prebench128: bench128.exe
	./bench128.exe pre

loadpmns:
	echo 'from commonpmns import pmnsdicts\nfrom commonpmns import primesdict\n\n(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdicts[$(PSIZE)128][primesdict[$(PSIZE)][$(INDEX)]]\nphi = 2**128' > pyparams128.py
	echo 'from commonpmns import pmnsdicts\nfrom commonpmns import primesdict\n\n(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdicts[$(PSIZE)][primesdict[$(PSIZE)][$(INDEX)]]\nphi = 2**64' > pyparams.py
	python3 precalcs128.py > params128.h $(PSIZE) $(INDEX)
	python3 precalcs.py > params.h $(PSIZE) $(INDEX)

loadpmnsWB:
	echo 'from commonpmns import pmnsWBdicts\nfrom commonpmns import primesdict\n\n(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsWBdicts[$(PSIZE)128][primesdict[$(PSIZE)][$(INDEX)]]\nphi = 2**128' > pyparams128.py
	echo 'from commonpmns import pmnsWBdicts\nfrom commonpmns import primesdict\n\n(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsWBdicts[$(PSIZE)][primesdict[$(PSIZE)][$(INDEX)]]\nphi = 2**64' > pyparams.py
	python3 precalcs128.py > params128.h
	python3 precalcs.py > params.h

progress: *.c *.py *.h makefile
	git add .
	git commit -m "progress"
	git push

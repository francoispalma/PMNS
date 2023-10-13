FLAGS= -Wall -Wextra -g -O3 -funswitch-loops -Wno-restrict -funroll-loops -fopenmp
CC = gcc-12
PSIZE = 1024
INDEX = 0
RDPMCFLAG=RDPMC_NOT_ALLOWED

all: main.exe

main.exe: main.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

hmain.exe: main.c hpmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns.o: pmns.c pmns.h params.h
	$(CC) -c $< $(FLAGS) 

params.h: precalcs.py pyparams.py 
	python3 $< > $@

structs.o: structs.c structs.h
	$(CC) -c $< $(FLAGS)

utilitymp.o: utilitymp.c utilitymp.h
	$(CC) -c $< $(FLAGS)

toeplitz.o: toeplitz.c
	$(CC) -c $< $(FLAGS)

toeplitz128.o: toeplitz128.c
	$(CC) -c $< $(FLAGS)

M192toep.o: M192toep.c
	$(CC) -c $< $(FLAGS)

p128.exe: main128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS)

pmns128.o: pmns128.c pmns128.h params128.h
	$(CC) -c $< $(FLAGS) -Wno-shift-negative-value

params128.h: precalcs128.py
	python3 $< > $@

pmns256.o: pmns256.c pmns256.h params256.h
	$(CC) -c $< $(FLAGS)

params256.h: precalcs256.py
	python3 $< > $@

hpmns.o: hpmns.c hpmns.h hparams.h
	$(CC) -c $< $(FLAGS)

hparams.h: hprecalcs.py
	python3 $< > $@

eccoptimizedcode.o: eccoptimizedcode.c eccoptimizedcode.h
	$(CC) -c $< $(FLAGS) 

gmpbench.exe: measuregmp.c
	$(CC) -o $@ $^ -O3 -Wall -g -lgmp -lcrypto -mavx512f -mavx512dq -mavx512vl -mavx512ifma -funroll-loops -funswitch-loops

bench.exe: intel-measurement.c pmns.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -D $(RDPMCFLAG)

bench128.exe: intel-measurement128.c pmns128.o structs.o utilitymp.o
	$(CC) -o $@ $^ $(FLAGS) -Wno-shift-negative-value

bench256.exe: pmns256.o
	$(CC) -o $@ $^ $(FLAGS)

hbench.exe: hbench.c hpmns.o structs.o utilitymp.o eccoptimizedcode.o
	$(CC) -o $@ $^ $(FLAGS) -lgmp

equalitytest.exe: equalitytest.c pmns.o structs.o utilitymp.o makefile
	$(CC) -o $@ equalitytest.c pmns.o structs.o utilitymp.o $(FLAGS) -ffast-math -ffinite-math-only -fno-signaling-nans -mfma -fno-tree-vectorize

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

hproof: hmain.exe
	./hmain.exe 100 > hlog
	python3 hproof.py

proof128: p128.exe
	./p128.exe 100 > log128
	python3 check128.py

p128: p128.exe
	./p128.exe 1

bench: bench.exe
	./bench.exe

hbench: hbench.exe
	./hbench.exe

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

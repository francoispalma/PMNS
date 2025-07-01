MAKEFLAGS += --no-print-directory
CC=gcc

main: mppmns.exe
	@make -B mppmns.exe 2>/dev/null
	./mppmns.exe
	python3 check.py

mppmns.exe: mppmns.c mpparams.h
	$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program

avx512mppmns.exe: avx512mppmns.c avx512params.h
	$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native

table5:
	@mv mpparams.h tmpparams.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/mpparams_1024_2.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_1024_3.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_1024_4.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@$(CC) -o gmpopenssl.exe gmpopenssl.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D LOG2P=1024 -lgmp -lcrypto && ./gmpopenssl.exe
	@echo
	@echo 2048 bit primes
	@cp tableparams/mpparams_2048_2.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_2048_3.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_2048_4.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@$(CC) -o gmpopenssl.exe gmpopenssl.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D LOG2P=2048 -lgmp -lcrypto && ./gmpopenssl.exe
	@echo
	@echo 4096 bit primes
	@cp tableparams/mpparams_4096_2.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_4096_3.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_4096_4.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@$(CC) -o gmpopenssl.exe gmpopenssl.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D LOG2P=4096 -lgmp -lcrypto && ./gmpopenssl.exe
	@echo
	@echo 6144 bit primes
	@cp tableparams/mpparams_6144_2.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_6144_3.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_6144_4.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@$(CC) -o gmpopenssl.exe gmpopenssl.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D LOG2P=6144 -lgmp -lcrypto && ./gmpopenssl.exe
	@echo
	@echo 8192 bit primes
	@cp tableparams/mpparams_8192_2.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_8192_3.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@cp tableparams/mpparams_8192_4.h mpparams.h
	@make -B mppmns.exe >/dev/null
	@./mppmns.exe
	@$(CC) -o gmpopenssl.exe gmpopenssl.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D LOG2P=8192 -lgmp -lcrypto && ./gmpopenssl.exe
	@rm mpparams.h
	@mv tmpparams.h mpparams.h 2>/dev/null || true

table6:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/avx512params_1024_2.h avx512params.h
	@echo s = 2
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_1024_3.h avx512params.h
	@echo s = 3
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 2048 bit primes
	@cp tableparams/avx512params_2048_2.h avx512params.h
	@echo s = 2
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_2048_3.h avx512params.h
	@echo s = 3
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 4096 bit primes
	@cp tableparams/avx512params_4096_2.h avx512params.h
	@echo s = 2
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_4096_3.h avx512params.h
	@echo s = 3
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_4096_4.h avx512params.h
	@echo s = 4
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 6144 bit primes
	@cp tableparams/avx512params_6144_2.h avx512params.h
	@echo s = 2
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_6144_3.h avx512params.h
	@echo s = 3
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_6144_4.h avx512params.h
	@echo s = 4
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 8192 bit primes
	@cp tableparams/avx512params_8192_2.h avx512params.h
	@echo s = 2
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_8192_3.h avx512params.h
	@echo s = 3
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@cp tableparams/avx512params_8192_4.h avx512params.h
	@echo s = 4
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

table8:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/avx512params_1024_3.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@$(CC) -o opensslavx512.exe opensslavx512.c rsaz-2k-avx512.s -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -march=native -D LOG2P=1024 -lgmp -lcrypto -flto && ./opensslavx512.exe
	@echo 1536 bit primes
	@cp tableparams/avx512params_1536_5.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@$(CC) -o opensslavx512.exe opensslavx512.c rsaz-3k-avx512.s -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -march=native -D LOG2P=1536 -lgmp -lcrypto -flto && ./opensslavx512.exe
	@echo 2048 bit primes
	@cp tableparams/avx512params_2048_6.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@$(CC) -o opensslavx512.exe opensslavx512.c rsaz-4k-avx512.s -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -march=native -D LOG2P=2048 -lgmp -lcrypto -flto && ./opensslavx512.exe
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

table9:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 807 bit primes
	@cp tableparams/avx512params_807_3.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 1214 bit primes
	@cp tableparams/avx512params_1214_2.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 1621 bit primes
	@cp tableparams/avx512params_1621_2.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 2029 bit primes
	@cp tableparams/avx512params_2029_2.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 2436 bit primes
	@cp tableparams/avx512params_2436_3.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 2844 bit primes
	@cp tableparams/avx512params_2844_3.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@echo 3251 bit primes
	@cp tableparams/avx512params_3251_3.h avx512params.h
	@make -B avx512mppmns.exe >/dev/null
	@./avx512mppmns.exe
	@echo
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

checktable5:
	@mv mpparams.h tmpparams.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/mpparams_1024_2.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=2 && python3 check.py
	@cp tableparams/mpparams_1024_3.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=3 && python3 check.py
	@cp tableparams/mpparams_1024_4.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=4 && python3 check.py
	@echo
	@echo 2048 bit primes
	@cp tableparams/mpparams_2048_2.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=2 && python3 check.py
	@cp tableparams/mpparams_2048_3.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=3 && python3 check.py
	@cp tableparams/mpparams_2048_4.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=4 && python3 check.py
	@echo
	@echo 4096 bit primes
	@cp tableparams/mpparams_4096_2.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=2 && python3 check.py
	@cp tableparams/mpparams_4096_3.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=3 && python3 check.py
	@cp tableparams/mpparams_4096_4.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=4 && python3 check.py
	@echo
	@echo 6144 bit primes
	@cp tableparams/mpparams_6144_2.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=2 && python3 check.py
	@cp tableparams/mpparams_6144_3.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=3 && python3 check.py
	@cp tableparams/mpparams_6144_4.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=4 && python3 check.py
	@echo
	@echo 8192 bit primes
	@cp tableparams/mpparams_8192_2.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=2 && python3 check.py
	@cp tableparams/mpparams_8192_3.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=3 && python3 check.py
	@cp tableparams/mpparams_8192_4.h mpparams.h
	@$(CC) -o mppmns.exe mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -D NOBENCH
	@./mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=4 && python3 check.py
	@echo
	@rm mpparams.h
	@mv tmpparams.h mpparams.h 2>/dev/null || true

checktable6:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/avx512params_1024_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=2 && python3 avx512check.py
	@cp tableparams/avx512params_1024_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=3 && python3 avx512check.py
	@echo
	@echo 2048 bit primes
	@cp tableparams/avx512params_2048_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=2 && python3 avx512check.py
	@cp tableparams/avx512params_2048_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=3 && python3 avx512check.py
	@echo
	@echo 4096 bit primes
	@cp tableparams/avx512params_4096_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=2 && python3 avx512check.py
	@cp tableparams/avx512params_4096_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=3 && python3 avx512check.py
	@cp tableparams/avx512params_4096_4.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=4096 && export NBCHUNKS=4 && python3 avx512check.py
	@echo
	@echo 6144 bit primes
	@cp tableparams/avx512params_6144_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=2 && python3 avx512check.py
	@cp tableparams/avx512params_6144_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=3 && python3 avx512check.py
	@cp tableparams/avx512params_6144_4.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=6144 && export NBCHUNKS=4 && python3 avx512check.py
	@echo
	@echo 8192 bit primes
	@cp tableparams/avx512params_8192_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=2 && python3 avx512check.py
	@cp tableparams/avx512params_8192_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=3 && python3 avx512check.py
	@cp tableparams/avx512params_8192_4.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=8192 && export NBCHUNKS=4 && python3 avx512check.py
	@echo
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

checktable8:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 1024 bit primes
	@cp tableparams/avx512params_1024_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1024 && export NBCHUNKS=3 && python3 avx512check.py
	@echo 1536 bit primes
	@cp tableparams/avx512params_1536_5.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1536 && python3 avx512check.py
	@echo 2048 bit primes
	@cp tableparams/avx512params_2048_6.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2048 && export NBCHUNKS=6 && python3 avx512check.py
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

checktable9:
	@mv avx512params.h tavx512params.h 2>/dev/null || true
	@echo 807 bit primes
	@cp tableparams/avx512params_807_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=807 && python3 avx512check.py
	@echo
	@echo 1214 bit primes
	@cp tableparams/avx512params_1214_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1214 && python3 avx512check.py
	@echo
	@echo 1621 bit primes
	@cp tableparams/avx512params_1621_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=1621 && python3 avx512check.py
	@echo
	@echo 2029 bit primes
	@cp tableparams/avx512params_2029_2.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2029 && python3 avx512check.py
	@echo
	@echo 2436 bit primes
	@cp tableparams/avx512params_2436_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2436 && python3 avx512check.py
	@echo
	@echo 2844 bit primes
	@cp tableparams/avx512params_2844_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=2844 && python3 avx512check.py
	@echo
	@echo 3251 bit primes
	@cp tableparams/avx512params_3251_3.h avx512params.h
	@$(CC) -o avx512mppmns.exe avx512mppmns.c -g -Wall -Wextra -O3 -funswitch-loops -funroll-loops -ftree-vectorize -fwhole-program -march=native -D NOBENCH
	@./avx512mppmns.exe
	@export PSIZE=3251 && python3 avx512check.py
	@echo
	@rm avx512params.h
	@mv tavx512params.h avx512params.h 2>/dev/null || true

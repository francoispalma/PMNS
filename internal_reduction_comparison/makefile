CC = gcc
FASTFLAGS = -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize
WARNINGS = -g -Wall -Wextra -Wno-comment
MAKEFLAGS += --no-print-directory
execf="BENCH"
exponf="PLANT"

redintcomp:
	@mv redparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_5.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 512 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_9.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_10.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 1024 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_18.h redparams.h
	@make -B redintcomp.exe execf="NOMONTG" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_20.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_21.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 2048 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_38.h redparams.h
	@make -B redintcomp.exe execf="NOMONTG" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_42.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_44.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 4096 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_80.h redparams.h
	@make -B redintcomp.exe execf="NOMONTG" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_84.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_100.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 6144 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_128.h redparams.h
	@make -B redintcomp.exe execf="NOMONTG" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_132.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_152.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 8192 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_168.h redparams.h
	@make -B redintcomp.exe execf="NOMONTG" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_184.h redparams.h
	@make -B redintcomp.exe execf="NOBARRETT" 2>/dev/null
	@./redintcomp.exe
	@cp tableparams/redparams_232.h redparams.h
	@make -B redintcomp.exe execf="BENCH" 2>/dev/null
	@./redintcomp.exe
	@rm redparams.h
	@mv tmpparams__.h redparams.h 2>/dev/null || true

checkredintcomp:
	@mv redparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes, n = 5
	@cp tableparams/redparams_5.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=5 && python3 redcheck.py
	@echo 512 bit primes, n = 9
	@cp tableparams/redparams_9.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=9 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 512 bit primes, n = 10
	@cp tableparams/redparams_10.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=10 && python3 redcheck.py
	@echo 1024 bit primes, n = 18
	@cp tableparams/redparams_18.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=18 && export NOMONTG=1 && python3 redcheck.py && unset NOMONTG
	@echo 1024 bit primes, n = 20
	@cp tableparams/redparams_20.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=20 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 1024 bit primes, n = 21
	@cp tableparams/redparams_21.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=21 && python3 redcheck.py
	@echo 2048 bit primes, n = 38
	@cp tableparams/redparams_38.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=38 && export NOMONTG=1 && python3 redcheck.py && unset NOMONTG
	@echo 2048 bit primes, n = 42
	@cp tableparams/redparams_42.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=42 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 2048 bit primes, n = 44
	@cp tableparams/redparams_44.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=44 && python3 redcheck.py
	@echo 4096 bit primes, n = 80
	@cp tableparams/redparams_80.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=80 && export NOMONTG=1 && python3 redcheck.py && unset NOMONTG
	@echo 4096 bit primes, n = 84
	@cp tableparams/redparams_84.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=84 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 4096 bit primes, n = 100
	@cp tableparams/redparams_100.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=100 && python3 redcheck.py
	@echo 6144 bit primes, n = 128
	@cp tableparams/redparams_128.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=128 && export NOMONTG=1 && python3 redcheck.py && unset NOMONTG
	@echo 6144 bit primes, n = 132
	@cp tableparams/redparams_132.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=132 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 6144 bit primes, n = 152
	@cp tableparams/redparams_152.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=152 && python3 redcheck.py
	@echo 8192 bit primes, n = 168
	@cp tableparams/redparams_168.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=168 && export NOMONTG=1 && python3 redcheck.py && unset NOMONTG
	@echo 8192 bit primes, n = 184
	@cp tableparams/redparams_184.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=184 && export NOBARRETT=1 && python3 redcheck.py && unset NOBARRETT
	@echo 8192 bit primes, n = 232
	@cp tableparams/redparams_232.h redparams.h
	@make -B redintcomp.exe execf="NOBENCH" 2>/dev/null
	@./redintcomp.exe
	@export NSIZE=232 && python3 redcheck.py
	@mv tmpparams__.h redparams.h 2>/dev/null || true

modexpcomp:
	@mv redparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_5.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 512 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_9.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_10.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 1024 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_18.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_20.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_21.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 2048 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_38.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_42.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_44.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 4096 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_80.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_84.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_100.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 6144 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_128.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_132.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_152.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@echo 8192 bit primes
	@echo -e '||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_168.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@echo -ne '||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_184.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="MONTG" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@cp tableparams/redparams_232.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" exponf="BARRETT" 2>/dev/null
	@echo -ne '\t||\t'
	@./redintcomp.exe
	@echo -e '\t||'
	@rm redparams.h
	@mv tmpparams__.h redparams.h 2>/dev/null || true


redintcomp.exe: redintcomp.c
	@$(CC) -o $@ $< $(WARNINGS) -D $(execf) -D $(exponf) $(FASTFLAGS) -fwhole-program -lgmp -fopenmp

int64vsint128.exe: int64vsint128.c
	@$(CC) -o $@ $< $(WARNINGS) $(FASTFLAGS)

check:
	python3 redcheck.py

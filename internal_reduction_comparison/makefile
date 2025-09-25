CC = gcc-12
MAKEFLAGS += --no-print-directory
execf="BENCH"

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

modexpcomp:
	@mv redparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_5.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 512 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_9.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 1024 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_18.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 2048 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_38.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 4096 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_80.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 6144 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_128.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@echo 8192 bit primes
	@echo -e '||\tn\t||\tPlant\t||\tMont\t||\tBarr\t||'
	@cp tableparams/redparams_168.h redparams.h
	@make -B redintcomp.exe execf="EXPBENCH" 2>/dev/null
	@./redintcomp.exe
	@rm redparams.h
	@mv tmpparams__.h redparams.h 2>/dev/null || true


redintcomp.exe: redintcomp.c
	@$(CC) -o redintcomp.exe redintcomp.c -D $(execf) -g -Wall -Wextra -Wno-comment -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize -fwhole-program -lgmp -fopenmp

check:
	python3 redcheck.py

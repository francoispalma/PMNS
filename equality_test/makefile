CC = gcc-12
FASTFLAGS = -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize
WARNINGS = -g -Wall -Wextra -Wno-comment
MAKEFLAGS += --no-print-directory
EQFLAG="INEQUALITY"

equalitytest.exe: equalitytest.c equparams.h
	@$(CC) -o $@ $< $(WARNINGS) $(FASTFLAGS) -D $(EQFLAG) -fwhole-program -lgmp

inequalitybench:
	@mv equparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes, N = 5
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_5.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 512 bit primes
	@echo N = 9
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_9.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 10
	@cp params/equparams_10.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 1024 bit primes
	@echo N = 20
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_20.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 21
	@cp params/equparams_21.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 2048 bit primes
	@echo N = 42
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_42.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 48
	@cp params/equparams_48.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 4096 bit primes
	@echo N = 84
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_84.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 99
	@cp params/equparams_99.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 6144 bit primes
	@echo N = 132
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_132.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 168
	@cp params/equparams_168.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 8192 bit primes
	@echo N = 184
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_184.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 232
	@cp params/equparams_232.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@rm equparams.h
	@mv tmpparams__.h equparams.h 2>/dev/null || true

equalitybench:
	@mv equparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes, N = 5
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_5.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 512 bit primes
	@echo N = 9
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_9.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 10
	@cp params/equparams_10.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 1024 bit primes
	@echo N = 20
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_20.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 21
	@cp params/equparams_21.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 2048 bit primes
	@echo N = 42
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_42.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 48
	@cp params/equparams_48.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 4096 bit primes
	@echo N = 84
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_84.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 99
	@cp params/equparams_99.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 6144 bit primes
	@echo N = 132
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_132.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 168
	@cp params/equparams_168.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 8192 bit primes
	@echo N = 184
	@echo -e '================================================================================'
	@echo -e '|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|\tHeuri\t|'
	@cp params/equparams_184.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo N = 232
	@cp params/equparams_232.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@rm equparams.h
	@mv tmpparams__.h equparams.h 2>/dev/null || true

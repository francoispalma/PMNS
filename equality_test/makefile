CC = gcc-12
FASTFLAGS = -O3 -funswitch-loops -funroll-loops -fno-tree-vectorize
WARNINGS = -g -Wall -Wextra -Wno-comment
MAKEFLAGS += --no-print-directory
EQFLAG="INEQUALITY"

equalitytest.exe: equalitytest.c equparams.h
	@$(CC) -o $@ $< $(WARNINGS) $(FASTFLAGS) -D $(EQFLAG) -fwhole-program -lgmp

hequalitytest.exe: linear_red/hequalitytest.c linear_red/hequparams.h
	@$(CC) -o $@ $< $(WARNINGS) $(FASTFLAGS) -D $(EQFLAG) -fwhole-program -lgmp

inequalitybench:
	@mv equparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_5.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 512 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_9.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_10.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 1024 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_20.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_21.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 2048 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_42.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_48.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 4096 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_84.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_99.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 6144 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_132.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_168.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 8192 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_184.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_232.h equparams.h
	@make -B equalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./equalitytest.exe
	@rm equparams.h
	@mv tmpparams__.h equparams.h 2>/dev/null || true

equalitybench:
	@mv equparams.h tmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_5.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 512 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_9.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_10.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 1024 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_20.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_21.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 2048 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_42.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_48.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 4096 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_84.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_99.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 6144 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_132.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_168.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@echo 8192 bit primes
	@echo -e '================================================================================'
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tReduc\t|\tExten\t|'
	@cp params/equparams_184.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@cp params/equparams_232.h equparams.h
	@make -B equalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./equalitytest.exe
	@rm equparams.h
	@mv tmpparams__.h equparams.h 2>/dev/null || true

hinequalitybench:
	@mv linear_red/hequparams.h linear_red/htmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_5.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 512 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_9.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 1024 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_18.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 2048 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_37.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 4096 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_74.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 6144 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_144.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 8192 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_180.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="INEQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@rm linear_red/hequparams.h
	@mv linear_red/htmpparams__.h linear_red/hequparams.h 2>/dev/null || true

hequalitybench:
	@mv linear_red/hequparams.h linear_red/htmpparams__.h 2>/dev/null || true
	@echo 256 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_5.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 512 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_9.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 1024 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_18.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 2048 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_37.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 4096 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_74.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 6144 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_144.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@echo 8192 bit primes
	@echo -e '================================================================='
	@echo -e '|\tN\t|\tNaive\t|\tTrans\t|\tCarry\t|'
	@cp linear_red/params/hequparams_180.h linear_red/hequparams.h
	@make -B hequalitytest.exe EQFLAG="EQUALITY" 2>/dev/null
	@./hequalitytest.exe
	@rm linear_red/hequparams.h
	@mv linear_red/htmpparams__.h linear_red/hequparams.h 2>/dev/null || true



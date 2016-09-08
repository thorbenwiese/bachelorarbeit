CC = gcc

CFLAGS = -std=c99 -Wall -Wextra -O2 -ggdb -lm

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

#Targets
all: tp_calc

test: tp_calc.c Makefile

	$(CC) $(CFLAGS) tp_calc.c -o tp_calc

clean:
	$(RM) *.o tp_calc
	$(RM) -r tp_calc.dSYM

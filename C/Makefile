CC?=gcc

CFLAGS=-Wall -Werror -O3 -g -Wunused-parameter -Wpointer-arith -fno-strict-aliasing -lm

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

#Targets
all: main

test: main.c Makefile

	$(CC) $(CFLAGS) main.c -o main

clean:
	$(RM) *.o main
	$(RM) -r tp_calc.dSYM

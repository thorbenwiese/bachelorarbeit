CC = gcc

CFLAGS = -std=c99 -Wall -Wextra -O2 -ggdb

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

#Targets
all: test

test: test.c Makefile

	$(CC) $(CFLAGS) test.c -o test

clean:
	$(RM) *.o test
	$(RM) -r test.dSYM

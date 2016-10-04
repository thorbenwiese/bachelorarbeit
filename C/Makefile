CC?=gcc
CFLAGS=-Wall -Werror -O3 -g -fno-strict-aliasing -Wpointer-arith

OBJ=front-with-trace.o tracer.o eoplist.o TracePoint.o

all:test.x

test:
	./test.x gagcatgttgcctggtcctttgctaggtactgtagaga gaccaagtaggcgtggaccttgctcggtctgtaagaga 0 37 4 37 10
	./test.x gagcatgttgcctggtcctttgctaggtactgtagaga gaccaagtaggcgtggaccttgctcggtctgtaagaga 2 23 0 22 5
	./test.x acagcatcatcagacgcacat acgagacctactcagcagac 5 20 4 19 5
	./test.x acagcatcatcagacgcacat acgagacctactcagcagac 1 20 3 20 5

# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

test.x:${OBJ} main.o
	${CC} main.o ${OBJ} -lm -o $@

# read dependencies
-include $(wildcard *.d)

.PHONY:clean
clean:
	$(RM) *.[oxd]

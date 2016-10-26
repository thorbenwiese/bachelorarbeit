CC?=gcc
CFLAGS=-Wall -Werror -O3 -g -fno-strict-aliasing -Wpointer-arith -Wunused-variable

OBJ=front-with-trace.o tracer.o eoplist.o TracePoint.o

all:test.x

test:test.x
	./test.x -f sequence_file.txt -p 0 400 0 400 -d 100 -a 1 -x
	./test.x -f sequence_file.txt -p 100 1000 100 1000 -d 200 -a 1 -x
	./test.x -f sequence_file.txt -p 666 777 650 777 -d 5 -a 1 -x
	./test.x -f sequence_file.txt -p 650 1000 650 1020 -d 50 -a 1 -x
	./test.x -f sequence_file.txt -p 12 99 13 100 -d 5 -a 1 -x

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

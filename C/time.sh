#!/bin/bash
for i in $(seq 5 15 500)
do
	./test.x -f sequence_file.txt -p $((RANDOM % 300+1)) 4999 $((RANDOM % 200+1)) 4999 -d $i -a 5000 -x
done

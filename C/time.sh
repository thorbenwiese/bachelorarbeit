#!/bin/bash
OUTPUT="$(ruby ./gen-randseq.rb -m pair -l 5000)"
for i in {5..100}
do
	for j in {1..10}
	do
		./test.x ${OUTPUT} $((RANDOM % 300+1)) 4999 $((RANDOM % 200+1)) 4999 $i
	done
done

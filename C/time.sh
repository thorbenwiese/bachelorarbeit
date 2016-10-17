#!/bin/bash
OUTPUT="$(ruby ./gen-randseq.rb -m pair -l 5000)"
for i in {5..10}
do
	for j in {1..50}
	do
		./test.x ${OUTPUT} 0 4999 0 4999 $i
	done
	echo $i
done

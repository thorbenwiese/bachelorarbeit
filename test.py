#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import alignment2
import neuer_versuch

""" Testfälle:

-v -rebuildonly --delta
random_sequences(amount, random_length, error_rate, alphabet)
random_sequences(100,50,0.15,"acgt")
random_sequences(10,200,0.30,"acgt")
random_sequences(30,100,0.7,"acgt")
random_sequences(10,50,0.15,"ac")

"""

def test_random_sequences(amount,random_length,error_rate,alphabet,delta):
	rseqs = neuer_versuch.random_sequences(amount,random_length,error_rate,alphabet)
	tpa = neuer_versuch.TracePointAlignment(rseqs[0],rseqs[1],delta)
	for i in range(0,amount*2,2):
            aln = tpa.calculate_alignment(rseqs[i],rseqs[i+1])
            #print aln.pretty_print_alignment(aln.seq1, aln.seq2)
            #print "Score:", aln.score
            #print ""
            score = aln.score
            alignment = neuer_versuch.TracePointAlignment(rseqs[i],rseqs[i+1],delta,score)
            alignment.calculate_intervals(aln)

    	return "OK"

def test_without_cigar(seq1,seq2,delta):
	aln = neuer_versuch.TracePointAlignment(seq1, seq2, delta)
	test = aln.calculate_alignment(seq1,seq2)
	print ""
	print "ALIGNMENT\n", aln.pretty_print_alignment(test.seq1,test.seq2)
	tp = aln.calculate_intervals(aln,True)
	print "TP",tp.tp

	return "OK"

def main():
	#tpa = neuer_versuch.TracePointAlignment(seq1,seq2,15)
	#print tpa.pretty_print_alignment(seq1,seq2)
	print "Testing random_sequences 1/4...", test_random_sequences(100,50,0.15,"acgt",5)
	print "Testing random_sequences 2/4...", test_random_sequences(10,200,0.30,"acgt",20)
	print "Testing random_sequences 3/4...", test_random_sequences(30,100,0.10,"acgt",15)
	print "Testing random_sequences 4/4...", test_random_sequences(10,50,0.15,"ac",5)
	print ""

	seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga"
	seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga"
	print "Testing sequences without CIGAR 1/3...", test_without_cigar(seq1,seq2,15)

	seq1 = "acgtgtggc"
	seq2 = "aaacgggcacgccgtggcccct"
	#print "Testing sequences without CIGAR 2/3...", test_without_cigar(seq1,seq2,10)

	seq1 = "acggacgttgacagtgtgacgtacgagacgtgtttgacagtgaccaagaatgttagag"
	seq2 = "aggctcggacgtacgagacgtgtttggctcgagagc"
	print "Testing sequences without CIGAR 3/3...", test_without_cigar(seq1,seq2,15)
	print ""

	"""
	seq1 = "acccccggtggct"
	seq2 = "acgcaccgtcgcg"
	cigar = "13M"

	print "Testing sequences with CIGAR 1/4...", 
	print "Testing sequences with CIGAR 1/4...", 
	print "Testing sequences with CIGAR 1/4...", 
	print "Testing sequences with CIGAR 1/4...", """
if __name__ == "__main__":
    main()

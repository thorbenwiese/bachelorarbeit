#!/usr/bin/env python
# -*- coding: utf-8 -*-

import alignment2
import neuer_versuch
import argparse

def test_random_sequences(amount,random_length,error_rate,alphabet,delta,verbose):

	rseqs = neuer_versuch.random_sequences(amount,random_length,error_rate,alphabet)
	tpa = neuer_versuch.TracePointAlignment(rseqs[0],rseqs[1],delta)
	for i in range(0,amount*2,2):
            aln = tpa.calculate_alignment(rseqs[i],rseqs[i+1])
            alignment = neuer_versuch.TracePointAlignment(rseqs[i],rseqs[i+1],delta,aln.score)
            alignment.calculate_intervals(aln,verbose)

    	return "OK"

def test_without_cigar(seq1,seq2,delta,verbose):

	aln = neuer_versuch.TracePointAlignment(seq1, seq2, delta)
	test = aln.calculate_alignment(seq1,seq2)
	tp = test.calculate_intervals(test,verbose)

	return "OK"

def test_with_cigar(seq1,seq2,cigar,delta,verbose):

	alignment = neuer_versuch.TracePointAlignment(seq1,seq2,delta)
	alignment1 = alignment.cigar_to_alignment(alignment.seq1, alignment.seq2, cigar)
	tp = alignment1.calculate_intervals(alignment1,verbose)

	return "OK"

def test_random(verbose):

	print "Testing random_sequences 1/4...", test_random_sequences(100,50,0.15,"acgt",5,verbose)
	print "Testing random_sequences 2/4...", test_random_sequences(10,200,0.30,"acgt",20,verbose)
	print "Testing random_sequences 3/4...", test_random_sequences(30,100,0.10,"acgt",15,verbose)
	print "Testing random_sequences 4/4...", test_random_sequences(10,50,0.15,"ac",5,verbose)
	print ""

def test_no_cigar(verbose):

	seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga"
	seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga"
	print "Testing sequences without CIGAR 1/3...", test_without_cigar(seq1,seq2,15,verbose)

	seq1 = "acgtgtggc"
	seq2 = "aaacgggcacgccgtggcccct"
	print "Testing sequences without CIGAR 2/3...", test_without_cigar(seq1,seq2,10,verbose)

	seq1 = "acggacgttgacagtgtgacgtacgagacgtgtttgacagtgaccaagaatgttagag"
	seq2 = "aggctcggacgtacgagacgtgtttggctcgagagc"
	print "Testing sequences without CIGAR 3/3...", test_without_cigar(seq1,seq2,15,verbose)
	print ""

def test_cigar(verbose):

	seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga" 
	seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga" 
	cigar = "5M5M1I3M5M1D9M1D5M1I2M2M"
	print "Testing sequences with CIGAR 1/4...", test_with_cigar(seq1,seq2,cigar,15,verbose)
	
	seq1 = "acccccggtggct"
	seq2 = "acgcaccgtcgcg"
	cigar = "13M"
	print "Testing sequences with CIGAR 2/4...", test_with_cigar(seq1,seq2,cigar,3,verbose)

	seq1 = "actgaactgact"
	seq2 = "actagaatggct"
	cigar = "3M1I3M1D5M"
	print "Testing sequences with CIGAR 3/4...", test_with_cigar(seq1,seq2,cigar,3,verbose)

	seq1 = "gtgtcgcccgtctagcatacgc"
	seq2 = "gggtgtaaccgactaggggg"
	cigar = "11M1D10M"
	print "Testing sequences with CIGAR 4/4...", test_with_cigar(seq1,seq2,cigar,3,verbose)

def test_bam(verbose):

	return True

def test_sam(verbose):

	return True

def test_rebuild(verbose):

	return True

def main():

	# TODO bisher nur Positiv-Tests
	# TODO rebuild_intervals einfügen

	parser = argparse.ArgumentParser()
	parser.add_argument("-a", "--all", help="Alle Tests durchführen", action="store_true")
	parser.add_argument("-n", "--nocigar", help="Teste Sequenzen ohne CIGAR-String", action="store_true")
	parser.add_argument("-c", "--cigar", help="Teste Sequenzen mit CIGAR-String", action="store_true")
	parser.add_argument("-b", "--bam", help="Test mit BAM-File", action="store_true")
	parser.add_argument("-s", "--sam", help="Test mit SAM-File", action="store_true")
	parser.add_argument("-r", "--random", help="Test mit zufälligen Sequenzen", action="store_true")
	parser.add_argument("-v", "--verbose", help="Ausführlicher Output", action="store_true")
	args = parser.parse_args()

	verbose = args.verbose

	if args.all:
		test_random(verbose)
		test_no_cigar(verbose)
		test_cigar(verbose)
		test_bam(verbose)
		test_sam(verbose)
		test_rebuild(verbose)
	elif args.nocigar:
		test_no_cigar(verbose)
	elif args.cigar:
		test_cigar(verbose)
	elif args.random:
		test_random(verbose)
	elif args.bam:
		test_bam(verbose)
	elif args.sam:
		test_sam(verbose)
	else:
		sys.stderr.write("Falsche Eingabe der Argumente!")
		sys.exit(1);


if __name__ == "__main__":
    main()

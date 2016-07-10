#!/usr/bin/env python

import tp_calc
import Alignment
import TracePoint

import argparse

# only for tests
start_seq1 = start_seq2 = 0

def test_random_sequences(amount,random_length,error_rate,alphabet,delta,
                          verbose,decode):

  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  for i in range(0, len(random_seq_list), 2):
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                    start_seq1, start_seq2)
    aln.calculate_alignment(random_seq_list[i], random_seq_list[i + 1])

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, 
                                            start_seq1, start_seq2)
    tp_aln.encode()

    if decode:
      tp_aln.decode(random_seq_list[i], random_seq_list[i + 1], delta, tp_aln.tp, 
                    start_seq1, start_seq2)

    if verbose:
      print "# TracePoints:", tp_aln.tp

def test_without_cigar(seq1,seq2,delta,verbose, decode):

  aln = Alignment.Alignment(seq1, seq2, start_seq1, start_seq2)
  aln.calculate_alignment(seq1, seq2)

  tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, 
                                          start_seq1, start_seq2)
  tp_aln.encode()

  if decode:
    tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1, start_seq2)

  if verbose:
    print "# TracePoints:", tp_aln.tp


def test_with_cigar(seq1,seq2,cigar,delta,verbose, decode):

  tp_aln = TracePoint.TracePointAlignment(seq1, seq2, delta, start_seq1, 
                                          start_seq2)
  tp_aln.encode_cigar(cigar)

  if verbose:
    print "# TracePoints:", tp_aln.tp

  if decode:
    tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1, start_seq2)

def test_random(verbose, decode):

  print "# Testing random_sequences 1/4...", test_random_sequences(80,50,0.15,"acgt",
                                                                   10,verbose, decode)
  print "# Testing random_sequences 2/4...", test_random_sequences(10,200,0.30,"acgt",
                                                                   20,verbose, decode)
  print "# Testing random_sequences 3/4...", test_random_sequences(30,100,0.10,"acgt",
                                                                   15,verbose, decode)
  print "# Testing random_sequences 4/4...", test_random_sequences(10,50,0.15,"ac",
                                                                   5,verbose, decode)
  print ""

def test_no_cigar(verbose, decode):
	
  seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga"
  seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga"
  print "# Testing sequences without CIGAR 1/3...", test_without_cigar(seq1,seq2,15,
                                                                       verbose, decode)
	
  seq1 = "acgtgtggc"
  seq2 = "aaacgggcacgccgtggcccct"
  print "# Testing sequences without CIGAR 2/3...", test_without_cigar(seq1,seq2,3,
                                                                       verbose,decode)
	
  seq1 = "acggacgttgacagtgtgacgtacgagacgtgtttgacagtgaccaagaatgttagag"
  seq2 = "aggctcggacgtacgagacgtgtttggctcgagagc"
  print "# Testing sequences without CIGAR 3/3...", test_without_cigar(seq1,seq2,15,
                                                                       verbose, decode)
  print ""
	
def test_cigar(verbose, decode):
	
  seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga" 
  seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga" 
  cigar = "5M5M1I3M5M1D9M1D5M1I2M2M"
  print "# Testing sequences with CIGAR 1/4...", test_with_cigar(seq1,seq2,cigar,10,
                                                                 verbose, decode)
	
  seq1 = "acccccggtggct"
  seq2 = "acgcaccgtcgcg"
  cigar = "13M"
  print "# Testing sequences with CIGAR 2/4...", test_with_cigar(seq1,seq2,cigar,3,
                                                                 verbose, decode)
	
  seq1 = "actgaactgact"
  seq2 = "actagaatggct"
  cigar = "3M1D3M1I5M"
  print "# Testing sequences with CIGAR 3/4...", test_with_cigar(seq1,seq2,cigar,3,
                                                                 verbose, decode)
	
  seq1 = "gtgtcgcccgtctagcatacgc"
  seq2 = "gggtgtaaccgactaggggg"
  cigar = "11M1D10M"
  print "# Testing sequences with CIGAR 4/4...", test_with_cigar(seq1,seq2,cigar,3,
                                                                 verbose, decode)
	
def test_bam(verbose, decode):
  print "TODO BAM"

def test_sam(verbose, decode):
  print "TODO SAM"

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--all", help="All tests", action="store_true")
  parser.add_argument("-n", "--nocigar", help="Test without CIGAR-String",
                      action="store_true")
  parser.add_argument("-c", "--cigar", help="Test with CIGAR-String",
                      action="store_true")
  parser.add_argument("-b", "--bam", help="Test BAM-File", action="store_true")
  parser.add_argument("-s", "--sam", help="Test SAM-File", action="store_true")
  parser.add_argument("-r", "--random", help="Test random sequences",
                      action="store_true")
  parser.add_argument("-v", "--verbose", help="Verbose output",default=False,
                      action="store_true")
  parser.add_argument("-d", "--decode",help="Test creating new alignment from TracePoints",
                      default=False, action="store_true" )

  args = parser.parse_args()

  verbose = args.verbose
  decode = args.decode

  if args.all:
    test_random(verbose, decode)
    test_no_cigar(verbose, decode)
    test_cigar(verbose, decode)
    # test_bam(verbose, decode)
    # test_sam(verbose, decode)
  elif args.nocigar:
    test_no_cigar(verbose, decode)
  elif args.cigar:
    test_cigar(verbose, decode)
  elif args.random:
    test_random(verbose, decode)
  elif args.bam:
    test_bam(verbose, decode)
  elif args.sam:
    test_sam(verbose, decode)
  else:
    sys.stderr.write("# Falsche Eingabe der Argumente!")
    sys.exit(1);

if __name__ == "__main__":
  main()

#!/usr/bin/env python

from random import randint
import tp_calc
import Alignment
import TracePoint
import Cigar_Pattern
import time

import argparse

# only for tests
start_seq1 = start_seq2 = 0

def test_random_sequences(amount,random_length,error_rate,alphabet,delta,
                          verbose,decode):

  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  for i in range(0, len(random_seq_list), 2):
    end_seq1 = len(random_seq_list[i])
    end_seq2 = len(random_seq_list[i+1]) 
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                    start_seq1, end_seq1, start_seq2, end_seq2)

    cigar = aln.calc_cigar(aln.seq1, aln.seq2)

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1, 
                                            end_seq1, start_seq2, end_seq2, 
                                            delta, cigar)

    if decode:
      cig = tp_aln.decode(tp_aln.tp)

      if verbose:
        aln1, aln2 = aln.cigar_to_aln(tp_aln.seq1, tp_aln.seq2, cig)
        aln.show_aln(aln1, aln2, cig)

    if verbose:
      print "# TracePoints:", tp_aln.tp

def test_without_cigar(seq1,seq2,delta,verbose, decode):

  end_seq1 = len(seq1)
  end_seq2 = len(seq2)
  aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1,
                            start_seq2, end_seq2)

  cigar = aln.calc_cigar(aln.seq1, aln.seq2)

  tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1, 
                                          end_seq1, start_seq2, end_seq2, 
                                          delta, cigar)

  if decode:
    cig = tp_aln.decode(tp_aln.tp)

    if verbose:
      aln1, aln2 = aln.cigar_to_aln(tp_aln.seq1, tp_aln.seq2, cig)
      aln.show_aln(aln1, aln2, cig)

  if verbose:
    print "# TracePoints:", tp_aln.tp


def test_with_cigar(seq1, seq2, cigar, delta, verbose, decode):

  end_seq1 = len(seq1)
  end_seq2 = len(seq2)
  tp_aln = TracePoint.TracePointAlignment(seq1, seq2,start_seq1, end_seq1,
                                          start_seq2,end_seq2, delta, cigar)

  aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1, 
                            start_seq2, end_seq2)

  if verbose:
    aln1, aln2 = aln.cigar_to_aln(tp_aln.seq1, tp_aln.seq2, cigar)
    aln.show_aln(aln1, aln2, cigar)
    print "# TracePoints:", tp_aln.tp

  if decode:
    cig = tp_aln.decode(tp_aln.tp)

    if verbose:
      aln1, aln2 = aln.cigar_to_aln(tp_aln.seq1, tp_aln.seq2, cig)
      aln.show_aln(aln1, aln2, cig)

def test_random(verbose, decode):

  print "# Testing random_sequences 1/4..." 
  test_random_sequences(80,50,0.15,"acgt",10,verbose,decode)
  
  print "# Testing random_sequences 2/4..."
  test_random_sequences(10,200,0.30,"acgt",20,verbose,decode)
  
  print "# Testing random_sequences 3/4..."
  test_random_sequences(30,100,0.10,"acgt",15,verbose,decode)
  
  print "# Testing random_sequences 4/4..."
  test_random_sequences(10,50,0.15,"ac",5,verbose,decode)
  print ""

def test_no_cigar(verbose, decode):
	
  seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga"
  seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga"
  
  print "# Testing sequences without CIGAR 1/3..."
  test_without_cigar(seq1,seq2,15,verbose,decode)
	
  seq1 = "acgtgtggc"
  seq2 = "aaacgggcacgccgtggcccct"
  
  print "# Testing sequences without CIGAR 2/3..."
  test_without_cigar(seq1,seq2,3,verbose,decode)
	
  seq1 = "acggacgttgacagtgtgacgtacgagacgtgtttgacagtgaccaagaatgttagag"
  seq2 = "aggctcggacgtacgagacgtgtttggctcgagagc"
  
  print "# Testing sequences without CIGAR 3/3..."
  test_without_cigar(seq1,seq2,15,verbose,decode)
  print ""
	
def test_cigar(verbose, decode):
	
  seq1 = "gagcatgttgcctggtcctttgctaggtactgtagaga" 
  seq2 = "gaccaagtaggcgtggaccttgctcggtctgtaagaga" 
  cigar = "10M1I8M1D9M1D5M1I4M"

  print "# Testing sequences with CIGAR 1/4..."
  test_with_cigar(seq1,seq2,cigar,10,verbose,decode)
	
  seq1 = "acccccggtggct"
  seq2 = "acgcaccgtcgcg"
  cigar = "13M"

  print "# Testing sequences with CIGAR 2/4..."
  test_with_cigar(seq1,seq2,cigar,3,verbose,decode)
	
  seq1 = "actgaactgact"
  seq2 = "actagaatggct"
  cigar = "3M1I3M1D5M"

  print "# Testing sequences with CIGAR 3/4..."
  test_with_cigar(seq1,seq2,cigar,3,verbose, decode)
	
  seq1 = "gtgtcgcccgtctagcatacgc"
  seq2 = "ggtgcgccgtcttagcata"
  cigar = "1M1D2M1I4M1D3M1I7M3D"

  print "# Testing sequences with CIGAR 4/4..."
  test_with_cigar(seq1,seq2,cigar,3,verbose,decode)

	
def main():

  randstring = ""
  for i in range(0,100):
    randstring += str(randint(0,9))

  print len(randstring)
  print randstring
  return



  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--all", help="All tests", action="store_true")
  parser.add_argument("-n", "--nocigar", help="Test without CIGAR-String",
                      action="store_true")
  parser.add_argument("-c", "--cigar", help="Test with CIGAR-String",
                      action="store_true")
  parser.add_argument("-r", "--random", help="Test random sequences",
                      action="store_true")
  parser.add_argument("-v", "--verbose", help="Verbose output",default=False,
                      action="store_true")
  parser.add_argument("-d", "--decode",
                      help="Test creating new alignment from TracePoints",
                      default=False, action="store_true")
  parser.add_argument("-i", "--intense",
                      help="Test with a lot of random sequences",
                      action="store_true")

  args = parser.parse_args()

  verbose = args.verbose
  decode = args.decode

  t = time.clock()

  if args.all:
    test_random(verbose, decode)
    test_no_cigar(verbose, decode)
    test_cigar(verbose, decode)
  elif args.nocigar:
    test_no_cigar(verbose, decode)
  elif args.cigar:
    test_cigar(verbose, decode)
  elif args.random:
    test_random(verbose, decode)
  elif args.intense:
    test_random_sequences(10000,200,0.15,"acgt",10,verbose,decode)
  else:
    sys.stderr.write("# Falsche Eingabe der Argumente!")
    sys.exit(1);

  print "Test complete.\nClock time: %.2f seconds." % (time.clock() - t)

if __name__ == "__main__":
  main()

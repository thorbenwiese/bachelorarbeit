#!/usr/bin/env python
# -*- coding: utf-8 -*-

import TracePoint
import Alignment

import sys
import math

# used for subprocesses in console
import subprocess
import argparse
import string, random
import os

# random sequence generator
def random_sequences(amount, random_length, error_rate, alphabet):
  seq = ""
  random_seqs = []
  for i in range(0, amount):
    seq = ''.join(random.choice(alphabet) for j in range(random_length))
    r_seq1 = list(seq)
    r_seq2 = list(seq)
        
    position = 0
    for element in r_seq2:
      # random number between 0 and 1
      r_num = random.random()
      if r_num <= error_rate:
        r_choice = random.random()
        if 0 <= r_choice < 0.1:
          # Insertion
          r_seq2.insert(position, '-')
        elif 0.1 < r_choice < 0.2:
          # Deletion
          r_seq1.insert(position, '-')
        elif 0.2 < r_choice < 0.4:
          # a
          if r_seq2[position] == 'a':
            continue
          r_seq2[position] = 'a'
        elif 0.4 < r_choice < 0.6:
          # c
          if r_seq2[position] == 'c':
            continue
          r_seq2[position] = 'c'
        elif 0.6 < r_choice < 0.8:
          # g
          if r_seq2[position] == 'g':
            continue
          r_seq2[position] = 'g'
        elif 0.8 < r_choice <= 1.0:
          # t
          if r_seq2[position] == 't':
            continue
          r_seq2[position] = 't'
        r_seq2[position] = random.choice(alphabet)
      position += 1
    seq_new1 = "".join(r_seq1)
    seq_new2 = "".join(r_seq2)
    random_seqs.append(seq_new1)
    random_seqs.append(seq_new2)

  return random_seqs

  # TODO
  def read_files(sequence_file, tp_file):
    print "TODO"

def main(argv):
  seq1 = seq2 = cigar = output = ""
  delta = 0

  # Arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-seq1", "--seq1", help="The first sequence")
  parser.add_argument("-start1", help="Starting position of the first sequence", type=int)
  parser.add_argument("-seq2", "--seq2", help="The second sequence")
  parser.add_argument("-start2", help="Starting position of the second sequence", type=int)
  parser.add_argument("-d", "--delta", help="Delta", type=int)
  group1 = parser.add_mutually_exclusive_group()
  group2 = parser.add_mutually_exclusive_group()
  group1.add_argument("-c", "--cigar", help="CIGAR-String")
  group1.add_argument("-b", "--bam", help="Input BAM-File")
  group1.add_argument("-s", "--sam", help="Input SAM-File")
  group1.add_argument("-r", "--random",
                        help="Random sequences generated with <Amount> <Length> <Error Rate> <Alphabet>",
                        nargs=4)
  group2.add_argument("-x", "--decode", help="Calculate alignment from TracePoint represantation", default=False,
                        action="store_true")
  
  args = parser.parse_args()

  seq1 = args.seq1
  start_seq1 = args.start1
  seq2 = args.seq2
  start_seq2 = args.start2
  cigar = args.cigar
  delta = args.delta

  # TODO to be embedded
  decode = args.decode

  # with CIGAR
  if args.cigar:
    tp_aln = TracePoint.TracePointAlignment(seq1, seq2, delta, start_seq1, start_seq2)
    tp_aln.encode_cigar(cigar)
    print "TP aus CIGAR:", tp_aln.tp

    if decode:
      print "WITH CIGAR\n",seq1,"\n",seq2
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1, start_seq2)

  # Random
  elif args.random:

    # only used for random sequences

    #start_seq1 = start_seq2 = 0

    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
    
    for i in range(0, len(random_seq_list), 2):
      aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], start_seq1, start_seq2)
      aln.calculate_alignment(random_seq_list[i], random_seq_list[i + 1])

      print "# Alignment:"
      print aln.show_aln(aln.seq1, aln.seq2)

      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, start_seq1, start_seq2)
      print "TEST1"
      tp_aln.encode()
      print "TEST2: ENCODE abgeschlossen!"

      if decode:
        # TODO replace ist nur für debugging, wird eh geändert
        print "RANDOM\n",random_seq_list[i],"\n",random_seq_list[i+1]
        tp_aln.decode(random_seq_list[i].replace("-",""), random_seq_list[i + 1].replace("-",""), delta, tp_aln.tp, start_seq1, start_seq2)

  # without CIGAR
  else:
    aln = Alignment.Alignment(seq1, seq2, start_seq1, start_seq2)
    aln.calculate_alignment(seq1, seq2)

    print "# Alignment:"
    print aln.show_aln(aln.seq1, aln.seq2)

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, start_seq1, start_seq2)
    tp_aln.encode()
   
    if decode:
      print "WITHOUT CIGAR\n",seq1,"\n",seq2
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1, start_seq2)

if __name__ == "__main__":
  main(sys.argv[1:])

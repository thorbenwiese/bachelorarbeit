#!/usr/bin/env python

import TracePoint
import Alignment

import sys
import math
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

    for i in range(0, len(r_seq1)):
      # random number between 0 and 1
      r_num = random.random()
      if r_num <= error_rate:
        r_choice = random.random()
        if 0 <= r_choice < 0.1:
          # Insertion
          r_seq2[position] = ''
        elif 0.1 < r_choice < 0.2:
          # Deletion
          r_seq1[position] = ''
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
      position += 1
    seq_new1 = "".join(r_seq1)
    seq_new2 = "".join(r_seq2)
    random_seqs.append(seq_new1)
    random_seqs.append(seq_new2)

  # store sequences
  for i in range(0,len(random_seqs)):
    
    if i == 0:
      with open('random_seq_file.txt', 'w') as file_:
        file_.write("%s\n" % random_seqs[i])
    else:
      with open('random_seq_file.txt', 'a') as file_:
        file_.write("%s\n" % random_seqs[i])

  # once input/output files are implemented, this should be deleted
  return random_seqs

# read sequences and alignment data from seperate files
def read_files(sequence_file, aln_file, id):

  sequences = []
  coded_aln = []

  with open(sequence_file, 'r') as seq_file:
    sequences = seq_file.readlines()

  with open(aln_file, 'r') as aln_file:
    coded_aln = aln_file.readlines()

  seq1 = seq2 = ""
  delta = start_seq1 = start_seq2 = 0
  tp = []

  seq1 = sequences[id*2-2]
  seq2 = sequences[id*2-1]

  delta, start_seq1, end_seq1, start_seq2, end_seq2, tp = coded_aln[id-1].split(";")

  aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1, start_seq2, end_seq2)
  aln_seq1, aln_seq2 = aln.calculate(aln.seq1, aln.seq2)
  aln.show_aln(aln_seq1.replace("\n",""), aln_seq2.replace("\n",""))

def main(argv):
  seq1 = seq2 = cigar = output = ""
  delta = 0

  parser = argparse.ArgumentParser()

  group1 = parser.add_mutually_exclusive_group()
  group2 = parser.add_mutually_exclusive_group()

  parser.add_argument("-seq1", "--seq1", help="The first sequence")
  parser.add_argument("-start1", help="Starting position of the first sequence", type=int)
  parser.add_argument("-end1", help="End position of the first sequence", type=int)
  parser.add_argument("-seq2", "--seq2", help="The second sequence")
  parser.add_argument("-start2", help="Starting position of the second sequence", type=int)
  parser.add_argument("-end2", help="End position of the second sequence", type=int)
  parser.add_argument("-d", "--delta", help="Delta", type=int)
  parser.add_argument("-iseq", "--input_seq", help="Input file with sequences")
  parser.add_argument("-ialn", "--input_aln", help="Input file with TracePoints")
  parser.add_argument("-id", "--id", help="Show specific alignment from Input files", type=int)

  group1.add_argument("-c", "--cigar", help="CIGAR-String")
  group1.add_argument("-b", "--bam", help="Input BAM-File")
  group1.add_argument("-s", "--sam", help="Input SAM-File")
  group1.add_argument("-r", "--random",
                        help="Random sequences generated with <Amount> <Length> <Error Rate> <Alphabet>",
                        nargs=4)

  # only while Input/Output files are not yet embedded
  group2.add_argument("-x", "--decode", help="Calculate alignment from TracePoint represantation", default=False,
                        action="store_true")
  
  args = parser.parse_args()


  start_seq1 = args.start1
  end_seq1 = args.end1
  start_seq2 = args.start2
  end_seq2 = args.end2
  seq1 = args.seq1
  seq2 = args.seq2
  cigar = args.cigar
  delta = args.delta

  decode = args.decode

  # show specific alignment from input files
  if args.input_seq and args.input_aln and args.id:
    read_files(args.input_seq, args.input_aln, args.id)
    sys.exit(1)

  if not args.random:
    seq1 = seq1[start_seq1:end_seq1]
    seq2 = seq2[start_seq2:end_seq2]

  # with CIGAR
  if args.cigar:

    tp_aln = TracePoint.TracePointAlignment(seq1, seq2, delta, cigar, start_seq1, end_seq1, 
                                            start_seq2, end_seq2)
    tp_aln.store_tp_aln('cigar_seq_file.txt','w')
    if decode:
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1,end_seq1,
                    start_seq2,end_seq2)

  # Random
  elif args.random:

    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])

    for i in range(0, len(random_seq_list), 2):

      aln = Alignment.Alignment(random_seq_list[i][start_seq1:end_seq1], random_seq_list[i + 1][start_seq2:end_seq2], 
                                start_seq1,end_seq1, start_seq2, end_seq2)

      aln_seq1, aln_seq2 = aln.calculate(aln.seq1, aln.seq2)
      cigar = aln.calc_cigar(aln_seq1, aln_seq2)
      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, cigar, start_seq1, 
                                              end_seq1, start_seq2, end_seq2)

      if i == 0:
        tp_aln.store_tp_aln('w')
      else:
        tp_aln.store_tp_aln('a')

      if decode:
        tp_aln.decode(random_seq_list[i], random_seq_list[i + 1], delta, tp_aln.tp, 
                      start_seq1, end_seq1, start_seq2, end_seq2)

  # without CIGAR
  else:
    aln = Alignment.Alignment(seq1, seq2, start_seq1,end_seq1, start_seq2,end_seq2)
    aln_seq1, aln_seq2 = aln.calculate(seq1, seq2)
    cigar = aln.calc_cigar(aln_seq1, aln_seq2)

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, cigar, start_seq1,
                                            end_seq1, start_seq2,end_seq2)
   
    if decode:
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1,end_seq1, 
                    start_seq2,end_seq2)

if __name__ == "__main__":
  main(sys.argv[1:])

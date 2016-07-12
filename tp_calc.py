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

    for i in range(0, len(r_seq2)):
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

  store_sequences(random_seqs)

  # once input/output files are implemented, this should be deleted
  return random_seqs

# store sequences to file
def store_sequences(seq_list):

  for i in range(0,len(seq_list)):

    with open('seq_file.txt', 'a') as file_:
      # without gaps from random sequences
      file_.write("%s\n" % seq_list[i].replace("-",""))

# read sequences and alignment data from seperate files
def read_files(sequence_file, aln_file):

  sequences = []
  coded_aln = []

  with open(sequence_file, 'r') as seq_file:
    sequences = seq_file.readlines()

  with open(aln_file, 'r') as aln_file:
    coded_aln = aln_file.readlines()

  return [sequences, coded_aln]

# encode and show specific alignment
def show_enc_aln(seq_file, aln_file, id):
  
  seq1 = seq2 = ""
  delta = start_seq1 = start_seq2 = 0
  tp = []

  input = read_files(seq_file, aln_file)

  seq1 = input[0][id*2-2]
  seq2 = input[0][id*2-1]

  data = input[1][id-1].split(";")

  delta = data[0]
  start_seq1 = data[1]
  start_seq2 = data[2]
  end_seq1 = data[3]
  end_seq2 = data[4]
  tp = data[5]

  tp_aln = TracePoint.TracePointAlignment(seq1, seq2, delta, start_seq1,end_seq1,start_seq2,end_seq2)
  tp_aln.tp = tp

  aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1, start_seq2, end_seq2)
  aln_seq1, aln_seq2 = aln.calculate(seq1, seq2)
  print aln.show_aln(aln_seq1.replace("\n",""), aln_seq2.replace("\n",""))

  # TODO TracePoints sind doch jetzt unwichtig?!
  print "# TracePoints:", tp_aln.tp

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
  parser.add_argument("-itp", "--input_tp", help="Input file with TracePoints")
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

  seq1 = args.seq1
  start_seq1 = args.start1
  end_seq1 = args.end1
  seq2 = args.seq2
  start_seq2 = args.start2
  end_seq2 = args.end2
  cigar = args.cigar
  delta = args.delta

  decode = args.decode

  # show specific alignment from input files
  if args.input_seq and args.input_tp and args.id:
    show_enc_aln(args.input_seq, args.input_tp, args.id)
    sys.exit(1)

  # with CIGAR
  if args.cigar:
    tp_aln = TracePoint.TracePointAlignment(seq1, seq2, delta, start_seq1, end_seq1, 
                                            start_seq2, end_seq1)
    tp_aln.encode_cigar(cigar)

    # test
    cig = tp_aln.calc_cigar(tp_aln.seq1, tp_aln.seq2)
    print tp_aln.seq1
    print tp_aln.seq2
    print "CIG 1:", cigar
    print "CIG 2:", cig
    if decode:
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1,end_seq1,
                    start_seq2,end_seq2)

  # Random
  elif args.random:

    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
    
    for i in range(0, len(random_seq_list), 2):
      aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], start_seq1, 
                                end_seq1, start_seq2, end_seq1)

      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, start_seq1, 
                                              end_seq1, start_seq2, end_seq2)
      aln_seq1, aln_seq2 = aln.calculate(random_seq_list[i], random_seq_list[i + 1])
      cig = tp_aln.calc_cigar(aln_seq1, aln_seq2)
      tp_aln.encode(cig)

      print "CIGAR:", cig
      print "TracePoints aus CIGAR:", tp_aln.tp
     
      tp_aln.store_tp_aln()

      if decode:
        tp_aln.decode(random_seq_list[i], random_seq_list[i + 1], delta, tp_aln.tp, 
                      start_seq1, end_seq1, start_seq2, end_seq2)

  # without CIGAR
  else:
    aln = Alignment.Alignment(seq1, seq2, start_seq1,end_seq1, start_seq2,end_seq2)
    aln.calculate_alignment(seq1, seq2)

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, delta, start_seq1,
                                            end_seq1, start_seq2,end_seq2)
    tp_aln.encode()
   
    if decode:
      tp_aln.decode(seq1, seq2, delta, tp_aln.tp, start_seq1,end_seq1, 
                    start_seq2,end_seq2)

if __name__ == "__main__":
  main(sys.argv[1:])

#!/usr/bin/env python

import TracePoint
import Alignment
import sys
import argparse
import string, random
import time


def mutate(sequence,err_prob,alphabet):
  seqlen = len(sequence)
  asize = len(alphabet)
  s = []
  i = 0
  for element in range(0,seqlen):
    r = random.random()
    if r <= err_prob:
      r = random.random()
      if r <= 0.8:
        s.append(random.choice(alphabet))
        i += 1
      elif r <= 0.9:
        s.append(random.choice(alphabet))
      else:
        i += 1
    else:
      s.append(sequence[i])
      i += 1
    if i >= seqlen:
     break
  return "".join(s)


# random sequence generator
def random_sequences(amount, random_length, error_rate, alphabet):

  assert (amount > 0), "Amount of random sequences should be > 0."
  assert (random_length > 0), "Length of random sequences should be > 0."
  assert (error_rate >= 0), "Error rate should not be negative."
  assert alphabet, "Alphabet must be set."

  random_seqs = []

  for i in range(0, amount):
    rseq1 = ''.join(random.choice(alphabet) for j in range(random_length))
    rseq2 = mutate(rseq1, error_rate, alphabet)
    random_seqs.append(rseq1)
    random_seqs.append(rseq2)

  assert random_seqs, "Random sequences are empty."

  # store sequences
  for i in range(0,len(random_seqs)):
    
    if i == 0:
      with open('random_seq_file.txt', 'w') as file_:
        file_.write("%s\n" % random_seqs[i])
    else:
      with open('random_seq_file.txt', 'a') as file_:
        file_.write("%s\n" % random_seqs[i])

  return random_seqs

# read sequences and alignment data from seperate files
def read_files(sequence_file, aln_file, id):

  sequences = []
  coded_aln = []

  with open(sequence_file, 'r') as seq_file:
    sequences = seq_file.readlines()

  with open(aln_file, 'r') as aln_file:
    coded_aln = aln_file.readlines()

  seq1 = sequences[id * 2 - 2].replace("\n","")
  seq2 = sequences[id * 2 - 1].replace("\n","")

  data = coded_aln[id - 1].split(";")
  delta = int(data[0])
  start_seq1 = int(data[1])
  end_seq1 = int(data[2])
  start_seq2 = int(data[3])
  end_seq2 = int(data[4])
  
  # tp is stored as a list converted to a string
  # it has to be converted to a list of ints again
  tp = map(int, data[5].replace("[","").replace("]","").split(','))
  
  tp_aln = TracePoint.TracePointAlignment(seq1, seq2, start_seq1, end_seq1, 
                                          start_seq2, end_seq2, delta)

  # reconstruct CIGAR-String from TracePoints and sequences
  cigar = tp_aln.decode(tp)
  aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1, 
                            start_seq2, end_seq2)
  aln.show_aln(seq1, seq2, cigar)
  
def main(argv):

  parser = argparse.ArgumentParser()

  group1 = parser.add_mutually_exclusive_group()
  group2 = parser.add_mutually_exclusive_group()

  parser.add_argument("-seq1", "--seq1", help="The 1st sequence")
  parser.add_argument("-start1", help="Starting position of the 1st sequence", 
                      type=int)
  parser.add_argument("-end1", help="End position of the 1st sequence", 
                      type=int)
  parser.add_argument("-seq2", "--seq2", help="The 2nd sequence")
  parser.add_argument("-start2", help="Starting position of the 2nd sequence", 
                      type=int)
  parser.add_argument("-end2", help="End position of the 2nd sequence", 
                      type=int)
  parser.add_argument("-d", "--delta", help="Delta value", type=int)
  parser.add_argument("-iseq", "--input_seq", help="Input file with sequences")
  parser.add_argument("-ialn", "--input_aln", 
                      help="Input file with TracePoints")
  parser.add_argument("-id", "--id", 
                      help="Show specific alignment from Input files", 
                      type=int)

  group1.add_argument("-c", "--cigar", help="CIGAR-String")
  group1.add_argument("-b", "--bam", help="Input BAM-File")
  group1.add_argument("-s", "--sam", help="Input SAM-File")
  group1.add_argument("-r", "--random",
    help="Random sequences with <Amount> <Length> <Error Rate> <Alphabet>",
    nargs=4)

  args = parser.parse_args()

  start_seq1 = args.start1
  end_seq1 = args.end1
  start_seq2 = args.start2
  end_seq2 = args.end2
  seq1 = args.seq1
  seq2 = args.seq2
  cigar = args.cigar
  delta = args.delta

  t = time.clock()

  # show specific alignment from input files
  if args.input_seq and args.input_aln and args.id:
    read_files(args.input_seq, args.input_aln, args.id)
    return

  if not args.random:
    seq1 = seq1[start_seq1:end_seq1]
    seq2 = seq2[start_seq2:end_seq2]

  # with CIGAR
  if args.cigar:

    tp_aln = TracePoint.TracePointAlignment(seq1, seq2, start_seq1, end_seq1, 
                                            start_seq2, end_seq2, delta, cigar)
    tp_aln.store_tp_aln('cigar_seq_file.txt','w')

  # Random
  elif args.random:
    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), 
                      float(args.random[2]), args.random[3])

    for i in range(0, len(random_seq_list), 2):

      aln = Alignment.Alignment(random_seq_list[i][start_seq1:end_seq1], 
                                random_seq_list[i + 1][start_seq2:end_seq2], 
                                start_seq1,end_seq1, start_seq2, end_seq2)

      cigar = aln.calc_cigar(aln.seq1, aln.seq2)
      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1,
                                              end_seq1, start_seq2, end_seq2, 
                                              delta, cigar)

      if i == 0:
        tp_aln.store_tp_aln('w')
      else:
        tp_aln.store_tp_aln('a')

  else: # args.cigar == false
    aln = Alignment.Alignment(seq1, seq2, start_seq1, end_seq1, 
                              start_seq2, end_seq2)
    cigar = aln.calc_cigar(seq1, seq2)

    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1, 
                                            end_seq1, start_seq2, end_seq2, 
                                            delta, cigar)

  print "Calculation complete.\nClock time: %.2f seconds." % (time.clock() - t)

 
if __name__ == "__main__":
  main(sys.argv[1:])

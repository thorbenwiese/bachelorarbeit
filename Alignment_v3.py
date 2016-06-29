#!/usr/bin/env python
# -*- coding: utf-8 -*-

import TracePoint_v3
import Cigar_v3

import sys
import math

# used for subprocesses in console
import subprocess
import argparse
import string, random
import os

# path for local submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

def calculate_alignment(sequence1, sequence2, delta):

  # x = No parameters.  Identical characters have score of 1, otherwise 0.
  # x = No gap penalties.
  alns = pairwise2.align.globalxx(sequence1.upper(), sequence2.upper())
  alignment = ""
  if len(alns) > 0:
    top_aln = alns[0]

    aln1, aln2, score, begin, end = top_aln

    alignment = TracePoint_v3.TracePointAlignment(aln1.lower(),aln2.lower(),delta,score)
  else:
    sys.stderr.write("# No alignment could be calculated.\n# One sequence is probably empty.\n")
    sys.exit(1)

  return alignment


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

# number of chars and '-'s in sequence  
def count_indels_letters(seq):

  letter_count = 0
  indel_count = 0
  for i in seq:
    if i.isalpha():
      letter_count += 1
    elif i == '-':
      indel_count += 1

  return [letter_count, indel_count]

# extract TracePoints from alignment
def alignment_to_tp(alignment,verbose=True):
  
  indels_in_seq2 = all_chars_in_seq1 = interval_count = 0

  # number of chars in sequences
  letters_in_seq1 = count_indels_letters(alignment.seq1)[0]
  letters_in_seq2 = count_indels_letters(alignment.seq2)[0]

  # Trace Point Array
  tp = []

  end_seq1 = alignment.start_seq1 + letters_in_seq1 - 1
  end_seq2 = alignment.start_seq2 + letters_in_seq2 - 1

  # dynamic calculation of interval size
  dynamic = max(1,int(math.ceil(alignment.start_seq1/alignment.delta)))

  # adjustment of interval length
  interval_count = min(int(math.ceil(float(len(alignment.seq1)) / alignment.delta)), int(math.ceil(float(len(alignment.seq2)) / alignment.delta)))

  # intervals
  # initialized with 0s
  intervals = [0] * interval_count

  for i in range(0, interval_count):
    # first interval
    if i == 0:
      intervals[i] = alignment.start_seq1, dynamic * alignment.delta - 1
    # last interval
    elif i == interval_count - 1:
      intervals[i] = (dynamic + interval_count - 2) * alignment.delta, end_seq1
    # other intervals
    else:
      intervals[i] = (dynamic + i - 1) * alignment.delta, (dynamic + i) * alignment.delta - 1


  start_seq1 = alignment.start_seq1
  start_seq2 = alignment.start_seq2
  start1 = start_seq1
  start2 = start_seq2
  for i in range(0,len(intervals)):
    
    # ignore last interval
    if i == len(intervals) - 2:
      break
    else:

      # Anzahl der Buchstaben an den Enden der Sequenzen
      rest1 = count_indels_letters(alignment.seq1[all_chars_in_seq1:end_seq1])[0]
      rest2 = count_indels_letters(alignment.seq2[all_chars_in_seq1:end_seq2])[0]

      # break condition
      if rest1 != 0 and rest2 != 0:
        # first delta letters in alignment.seq1
        while count_indels_letters(alignment.seq1[start_seq1:all_chars_in_seq1])[0] != alignment.delta:
          all_chars_in_seq1 += 1

      # count gaps in seq2
      indel_count = count_indels_letters(alignment.seq2[start_seq1:all_chars_in_seq1])[1]
      indels_in_seq2 += indel_count
      # calculate Trace Point
      tp.append(all_chars_in_seq1 - indels_in_seq2 - 1)
      end = all_chars_in_seq1 - indels_in_seq2 - 1
      start_seq1 = all_chars_in_seq1
      all_chars_in_seq1 = start_seq1
      start1 = intervals[i][1] + 1
      start2 = start_seq1 - indels_in_seq2
      end = all_chars_in_seq1 - indels_in_seq2 - 1
  
  TP_alignment = TracePoint_v3.TracePointAlignment(alignment.seq1, alignment.seq2, alignment.delta, tp=tp)

  if verbose:
    print "# Trace Points:", tp, "\n"
  
  return TP_alignment

# 0s and 5s for the pretty_print_alignment
def print_sequence_positions(seq):

  positions = ""
  count = False
  for i in range(0, len(seq), 5):
    if not count:
      positions += "0" + " " * 4
      count = True
    else:
      positions += "5" + " " * 4
      count = False

  return positions


# pretty_print of an alignment with positions and '|'s
def show_aln(seq1, seq2):

  middle = pretty_print= ""
  count = add = 0

  if len(seq1) < len(seq2):
    add = len(seq2) - len(seq1)
  else:
    add = len(seq1) - len(seq2)

  seq2 += " " * add
  for i in seq1:
    if i == seq2[count]:
      middle += '|'
    else:
      middle += ' '
    count += 1

  pos1 = print_sequence_positions(seq1)
  pos2 = print_sequence_positions(seq2)
  pretty_print = "\n" + pos1 + "\n" + seq1 + "\n" + middle + "\n" + seq2 + "\n" + pos2 + "\n"

  return pretty_print

def main(argv):
  seq1 = seq2 = cigar = output = ""
  delta = 0

  # Arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-seq1", "--seq1", help="The first sequence")
  parser.add_argument("-start_seq1", help="Starting position of the first sequence", type=int, default=0)
  parser.add_argument("-seq2", "--seq2", help="The second sequence")
  parser.add_argument("-start_seq2", help="Starting position of the second sequence", type=int, default=0)
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
  parser.add_argument("-v", "--verbose", help="Verbose output", default=False, action="store_true")
  args = parser.parse_args()

  seq1 = args.seq1
  start_seq1 = args.start_seq1
  seq2 = args.seq2
  start_seq2 = args.start_seq2
  cigar = args.cigar
  delta = args.delta
  verbose = args.verbose

  # with CIGAR
  if args.cigar:
    cig_aln = Cigar_v3.CigarAlignment(seq1, seq2, delta, cigar, start_seq1, start_seq2)
    tp_aln = cig_aln.cigar_to_tp(verbose)

    if decode:
      tp_aln.tp_to_alignment(seq1,seq2,delta,tp_aln.tp)

  # Random
  elif args.random:
    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
    for i in range(0, len(random_seq_list), 2):
      aln = calculate_alignment(random_seq_list[i], random_seq_list[i + 1], delta)
      if verbose:
        print "# Alignment:"
        print show_aln(aln.seq1, aln.seq2)
      tp_aln = alignment_to_tp(aln)

      tp_aln.tp_to_alignment(random_seq_list[i], random_seq_list[i + 1], delta, tp_aln.tp)

  # without CIGAR
  else:
    aln = calculate_alignment(seq1, seq2, delta)
    if verbose:
      print "# Alignment:"
      print show_aln(aln.seq1, aln.seq2)
    tp_aln = alignment_to_tp(aln)
   
    # new alignment from TracePoints
    tp_aln.tp_to_alignment(seq1, seq2, delta, tp_aln.tp)

if __name__ == "__main__":
  main(sys.argv[1:])

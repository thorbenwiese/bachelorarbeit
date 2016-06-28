#!/usr/bin/env python
# -*- coding: utf-8 -*-

import TracePoint_v2
import Cigar_v2

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

    alignment = TracePoint_v2.TracePointAlignment(aln1.lower(),aln2.lower(),delta,score)
    # alignment = aln1 + "\n" + aln2 + "\n"
  else:
    sys.stderr.write("# No alignment was calculated.\n")
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
def alignment_to_tp(alignment):
  
  v_id = count = interval_count = 0

  # number of chars in sequences
  letters_in_seq1 = count_indels_letters(alignment.seq1)[0]
  letters_in_seq2 = count_indels_letters(alignment.seq2)[0]

  tp = []
  seq1_new = seq1_new = ""  

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


  seq1_new = seq2_new = ""
  start_seq1 = alignment.start_seq1
  start_seq2 = alignment.start_seq2
  start1 = start_seq1
  start2 = start_seq2
  for i in range(0,len(intervals)):
    
    if i == len(intervals) - 2:
      print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, letters_in_seq1 - 1, start2, letters_in_seq2 - 1)
      print show_aln(alignment.seq1[start_seq1:len(alignment.seq1)], alignment.seq2[start_seq1:len(alignment.seq1)])
      seq1_new += str(alignment.seq1[start_seq1:len(alignment.seq1)])
      seq2_new += str(alignment.seq2[start_seq1:len(alignment.seq2)])
    
    # if i != len(intervals) - 1:
    elif i == len(intervals) - 1:
      break
    else:

      # Anzahl der Buchstaben an den Enden der Sequenzen
      rest1 = count_indels_letters(alignment.seq1[count:end_seq1])[0]
      rest2 = count_indels_letters(alignment.seq2[count:end_seq2])[0]

      # Abbruchbedingung
      if rest1 != 0 and rest2 != 0:
        # first delta letters in alignment.seq1
        while count_indels_letters(alignment.seq1[start_seq1:count])[0] != alignment.delta:
          count += 1

      indel_count = 0
      indel_count = count_indels_letters(alignment.seq2[start_seq1:count])[1]
      v_id += indel_count
      tp.append(count - 1 - v_id)
      end = count - 1 - v_id
      indel_count = 0
      print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, intervals[i][1], start2, end)
      print show_aln(alignment.seq1[start_seq1:count], alignment.seq2[start_seq1:count])
      print ""
      seq1_new += str(alignment.seq1[start_seq1:count])
      seq2_new += str(alignment.seq2[start_seq1:count])
      start_seq1 = count
      count = start_seq1
      start1 = intervals[i][1] + 1
      start2 = start_seq1 - v_id
      end = count - 1 - v_id
  
  TP_alignment = TracePoint_v2.TracePointAlignment(seq1_new, seq2_new, alignment.delta, tp=tp)

  # TODO improve
  """
  start = end = 0
  c = 0

  for i in intervals:

    # number of remaining chars in sequences
    rest1 = count_indels_letters(alignment.seq1[end:end_seq1])[0]
    rest2 = count_indels_letters(alignment.seq2[end:end_seq2])[0]
    print "REST", rest1, rest2 
    # all but the last interval with break condition
    if i != intervals[-1] and rest1 != 0 and rest2 != 0:
      c+=1

      end = count_indels_letters(alignment.seq1[start:c * alignment.delta])[0]+count_indels_letters(alignment.seq1[start:c * alignment.delta])[1]
      print alignment.seq1[start:end]
      print alignment.seq2[start:end]
      print ""
      start = end + 1
      print "START", start, "END", end
      
      while alignment.delta != count_indels_letters(alignment.seq1[start:end])[0] and end != 30:
        print count_indels_letters(alignment.seq1[start:end])[0]
        end += 1
      tp.append(end)
      print "TPP", tp
      start = end + 1
      
  """     
  print "TracePoints:", tp
  
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

# calculation with random sequence generator
def random_calc(random_seq_list, delta, verbose, decode):
  for i in range(0, len(random_seq_list), 2):
    tp_alignment = TracePoint_v2.TracePointAlignment(random_seq_list[i], random_seq_list[i + 1], delta)
    if verbose:
      print "\n# Alignment of random sequences:"
      print tp_alignment.show_aln(tp_alignment.seq1, tp_alignment.seq2)
      #print "# Score:", tp_alignment.score
      print ""
    tp = tp_alignment.create_tp_aln(tp_alignment, verbose)

    # Output to file

    output = []
    output.extend((tp.seq1, tp.seq2, delta, tp.start_seq1, tp.start_seq2, (tp.tp)))

    if i == 0:
      # create and write new file
      tp.store_aln(output, 'w')
    elif i == -1:
      # append last output and rebuild if 'decode' flag is set
      tp.store_aln(output, 'a')
      if decode:
        tp.rebuild_intervals(tp, verbose)

      # TODO noch nicht sinnvoll eingebunden
      # read from file
      tp.read_and_create_aln("alignment_compressed.txt", verbose)
    else:
      # append other output to existing file
      tp.store_aln(output, 'a')


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

  if args.cigar:
    cig_aln = Cigar_v2.CigarAlignment(seq1, seq2, delta, cigar, start_seq1, start_seq2)
    cig_aln.cigar_to_tp()
  elif args.random:
    random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
    for i in range(0, len(random_seq_list), 2):
      tp_alignment = TracePoint_v2.TracePointAlignment(random_seq_list[i], random_seq_list[i + 1], delta)
      print show_aln(tp_alignment.seq1, tp_alignment.seq2)
      tp_aln = alignment_to_tp(tp_alignment)

      tp_aln.tp_to_alignment(random_seq_list[i], random_seq_list[i + 1], delta, tp_aln.tp)

  else:
    aln = calculate_alignment(seq1, seq2, delta)
    print show_aln(aln.seq1, aln.seq2)
    tp_aln = alignment_to_tp(aln)
   
    # new alignment from TracePoints
    tp_aln.tp_to_alignment(seq1, seq2, delta, tp_aln.tp)

if __name__ == "__main__":
  main(sys.argv[1:])

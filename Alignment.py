# -*- coding: utf-8 -*-

import sys

# path for local submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

# Alignment structure with CIGAR-String
class Alignment(object):

  # constructor
  def __init__(self, seq1, seq2,  start_seq1, start_seq2):
    self.seq1 = seq1
    self.seq2 = seq2
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2

  def calculate_alignment(sequence1, sequence2):

    # x = No parameters.  Identical characters have score of 1, otherwise 0.
    # x = No gap penalties.
    alns = pairwise2.align.globalxx(sequence1.upper(), sequence2.upper())
    alignment = ""
    if len(alns) > 0:
      top_aln = alns[0]

      aln1, aln2, score, begin, end = top_aln

     # kein TPAln..
     # alignment = TracePoint_v3.TracePointAlignment(aln1.lower(),aln2.lower())
    else:
      sys.stderr.write("# No alignment could be calculated.\n# One sequence is probably empty.\n")
      sys.exit(1)

    return alignment

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


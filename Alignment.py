import sys

# path for local submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Alignment structure with CIGAR-String
class Alignment(object):

  # constructor
  def __init__(self, seq1, seq2, start_seq1, start_seq2):
    self.seq1 = seq1
    self.seq2 = seq2
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2

  def calculate_alignment(self, sequence1, sequence2):

    # x = No parameters.  Identical characters have score of 1, otherwise 0.
    # x = No gap penalties.
    alns = pairwise2.align.globalxx(sequence1.upper(), sequence2.upper())
    alignment = ""
    if len(alns) > 0:
      top_aln = alns[0]

      aln1, aln2, score, begin, end = top_aln

      self.seq1 = aln1.lower()
      self.seq2 = aln2.lower()

    else:
      sys.stderr.write("# No alignment could be calculated.\n")
      sys.exit(1)


  # 0s and 5s for the pretty_print_alignment
  def print_sequence_positions(self, seq):

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
  def show_aln(self, seq1, seq2):

    middle = pretty_print= ""
    extend = 0

    if len(seq1) < len(seq2):
      extend = len(seq2) - len(seq1)
    else:
      extend = len(seq1) - len(seq2)

    seq2 += " " * extend

    for i in range(0,len(seq1)):
      if seq1[i] == seq2[i]:
        middle += '|'
      else:
        middle += ' '

    pos1 = self.print_sequence_positions(seq1)
    pos2 = self.print_sequence_positions(seq2)
    pretty_print = "\n" + pos1 + "\n" + seq1 + "\n" + middle + "\n" + seq2 + "\n" + pos2 + "\n"

    return pretty_print


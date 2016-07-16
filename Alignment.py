import sys

# path for local submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2

# Alignment structure with CIGAR-String
class Alignment(object):

  # constructor
  def __init__(self, seq1, seq2, start_seq1, end_seq1, start_seq2, end_seq2):
    self.seq1 = seq1
    self.seq2 = seq2
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2
    self.end_seq1 = end_seq1
    self.end_seq2 = end_seq2

  # calculate alignment with BioPython
  def calculate(self, sequence1, sequence2):

    assert sequence1, "First sequence for calculating the alignment is empty!"
    assert sequence2, "Second sequence for calculating the alignment is empty!"

    # m = A match score is the score of identical chars, otherwise mismatch score.
    # s = Same open and extend gap penalties for both sequences.
    alns = pairwise2.align.globalms(sequence1.upper(), sequence2.upper(),
                                    2, -1, -2, -1)
    assert (len(alns) > 0), "No alignment could be calculated.\n"
    top_aln = alns[0]

    aln1, aln2, score, begin, end = top_aln
    
    assert (len(aln1) == len(aln2)), "Alignment sequences do not have the same size!"

    return [aln1.lower(), aln2.lower()]

  # calculate CIGAR-String from two aln_seq
  def calc_cigar(self, aln1, aln2):

    assert aln1, "First alignment sequence for calculating CIGAR-String is empty."
    assert aln2, "Second alignment sequence for calculating CIGAR-String is empty."

    cigar = ""
    count = 1
    prev_match = False
    prev_ins = False
    prev_dele = False

    if aln1[0] == '-':
      previous_op  = 'D'
    elif aln2[0] == '-':
      previous_op  = 'I'
    else:
      previous_op = 'M'

    for i in range(1, len(aln1)):

      # deletion
      if aln1[i] == '-':
        if previous_op != 'D':
          cigar += "%d%s" % (count, previous_op)
          count = 1
          previous_op = 'D'
        else:
          count += 1

      # insertion
      elif aln2[i] == '-':
        if previous_op != 'I':
          cigar += "%d%s" % (count, previous_op)
          count = 1
          previous_op = 'I'
        else:
          count += 1

      # match / mismatch
      else:
        if previous_op != 'M':
          cigar += "%d%s" % (count, previous_op)
          count = 1
          previous_op = 'M'
        else:
          count += 1

    # last operation
    cigar += "%d%s" % (count, previous_op)

    return cigar

  # 0s and 5s for the pretty_print_alignment
  def print_sequence_positions(self, seq):

    assert seq, "Sequence for print_sequence_positions is empty."

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
  def show_aln(self, aln1, aln2):

    assert aln1, "First sequence for show_aln is empty."
    assert aln2, "Second sequence for show_aln is empty."

    middle = pretty_print= ""
    extend = 0

    for i in range(0,len(aln1)):
      middle += '|' if aln1[i] == aln2[i] else ' '

    print self.print_sequence_positions(aln1), "\n", aln1, "\n", middle, "\n", \
          aln2, "\n", self.print_sequence_positions(aln2)

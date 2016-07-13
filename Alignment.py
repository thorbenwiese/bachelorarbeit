import sys

# path for local submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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

    # x = No parameters.  Identical characters have score of 1, otherwise 0.
    # x = No gap penalties.
    alns = pairwise2.align.globalxx(sequence1.upper(), sequence2.upper())
    if len(alns) > 0:
      top_aln = alns[0]

      aln1, aln2, score, begin, end = top_aln

      return [aln1.lower(), aln2.lower()]

    else:
      sys.stderr.write("# No alignment could be calculated.\n")
      sys.exit(1)

  # calculate CIGAR-String from two aln_seq
  def calc_cigar(self, aln_seq1, aln_seq2):

    if len(aln_seq1) < len(aln_seq2):
      seqlen = len(aln_seq1)
    else:
      seqlen = len(aln_seq2)

    cigar = ""
    count = 0
    match = False
    ins = False
    dele = False

    for i in range(0, seqlen):

      # match
      if aln_seq1[i] == aln_seq2[i]:
        match = True
        if ins:
          ins = False
          cigar += "%d%s" % (count, 'I')
          count = 1
        elif dele:
          dele = False
          cigar += "%d%s" % (count, 'D')
          count = 1
        else:
          count += 1

      # deletion
      elif aln_seq1[i] == '-':

        dele = True
        if ins:
          ins = False
          cigar += "%d%s" % (count, 'I')
          count = 1
        elif match:
          match = False
          cigar += "%d%s" % (count, 'M')
          count = 1
        else:
          count += 1

      # insertion
      elif aln_seq2[i] == '-':

        ins = True
        if match:
          match = False
          cigar += "%d%s" % (count, 'M')
          count = 1
        elif dele:
          dele = False
          cigar += "%d%s" % (count, 'D')
          count = 1
        else:
          count += 1

    # last operation
    if match:
      cigar += "%d%s" % (count, 'M')
    elif dele:
      cigar += "%d%s" % (count, 'D')
    elif ins:
      cigar += "%d%s" % (count, 'I')

    return cigar

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

    print self.print_sequence_positions(seq1), "\n", seq1, "\n", middle, "\n", \
          seq2, "\n", self.print_sequence_positions(seq2)


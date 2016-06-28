# -*- coding: utf-8 -*-

import Alignment_v3

import re
import math
import sys

class TracePointAlignment(object):

  def __init__(self, seq1, seq2, delta, score=None, tp=None, start_seq1=0, start_seq2=0):
    self.seq1 = seq1
    self.seq2 = seq2
    self.delta = delta
    self.score = score
    self.tp = tp
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2

  # create new intervals from TracePoints and calculate new alignment
  def tp_to_alignment(self, seq1, seq2, delta, tp, start_seq1=0, start_seq2=0, verbose=True):
    count = 0
    new_seq1 = new_seq2 = ""
    for i in tp:
      if i == tp[0]:
        count += 1
        aln = Alignment_v3.calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
      elif i == tp[-1]:
        aln = Alignment_v3.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        count += 1
        aln = Alignment_v3.calculate_alignment(seq1[count*delta:len(seq1)],seq2[tp[count-1]+1:len(seq2)],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        break
      else:
        aln = Alignment_v3.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        count += 1
    if verbose:
      print "# Konkateniertes Alignment:"
      print Alignment_v3.show_aln(new_seq1, new_seq2)
    return True
  
  # store TracePointAlignment to file
    def store_aln(self, output, mode):
  
      with open('alignment_compressed.txt', mode) as file_:
  
        # ';' for easy splitting
        for item in output:
          file_.write("%s;" % item)

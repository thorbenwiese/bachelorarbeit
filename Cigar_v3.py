# -*- coding: utf-8 -*-
import math
import re
import sys

import TracePoint_v3

# Alignment structure with CIGAR-String
class CigarAlignment(object):

  # constructor
  def __init__(self, seq1, seq2, delta, cigar, start_seq1 = 0, start_seq2 = 0):
    self.seq1 = seq1
    self.seq2 = seq2
    self.delta = delta
    self.cigar = cigar
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2

  # extract TracePoints from CIGAR-String
  def cigar_to_tp(self,verbose):

    count = count1 = count2 = cig_count = 0

    end_seq1 = self.start_seq1 + len(self.seq1) - 1
    end_seq2 = self.start_seq2 + len(self.seq2) - 1

    interval_count = 0
    
    # dynamic calculation of interval size
    dynamic = max(1,int(math.ceil(self.start_seq1/self.delta)))

    # adjustment of interval length
    interval_count = min(int(math.ceil(float(len(self.seq1)) / self.delta)), int(math.ceil(float(len(self.seq2)) / self.delta)))
  
    # intervals 
    # initialized with 0s
    intervals = [0] * interval_count

    for i in range(0, interval_count):
      # first interval
      if i == 0:
        intervals[i] = self.start_seq1, dynamic * self.delta - 1
      # last interval
      elif i == interval_count - 1:
        intervals[i] = (dynamic + interval_count - 2) * self.delta, end_seq1
      # other intervals
      else:
        intervals[i] = (dynamic + i - 1) * self.delta, (dynamic + i) * self.delta - 1

    # create pattern for CIGAR-String
    cigar_pattern = re.compile(r"\d+[MIDNSHP=j]{1}")

    # search cigar for pattern
    tp = []
    for j in cigar_pattern.findall(self.cigar):
      cig_count = int(j[:-1])
      cig_symbol = j[-1]
    
      for i in range(0,cig_count):
        if cig_symbol == 'I':
          count1 += 1
        elif cig_symbol == 'D':
          count2 += 1
        else:
          count1 += 1
          count2 += 1
        # count until the end but ignore end of last interval as TracePoint
        if count1 == intervals[count][1] + 1 and count1 != len(self.seq1):
          tp.append(count2 - 1)
          count += 1
   
    if verbose:
      print "Trace Points:", tp

    tp_aln = TracePoint_v3.TracePointAlignment(self.seq1, self.seq2, self.delta, tp=tp, start_seq1=self.start_seq1, start_seq2=self.start_seq2)
    return tp_aln

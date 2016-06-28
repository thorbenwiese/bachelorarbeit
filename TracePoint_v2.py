# -*- coding: utf-8 -*-

import Alignment_v2

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
  def tp_to_alignment(self, seq1, seq2, delta, tp, start_seq1=0, start_seq2=0):
    count = 0
    check_aln1 = check_aln2 = ""
    for i in tp:
      if i == tp[0]:
        print start_seq1, "...", delta - 1
        print start_seq2, "...", i
        count += 1
        aln = Alignment_v2.calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i+1],delta)
        print Alignment_v2.show_aln(aln.seq1, aln.seq2)
        # print Alignment_v2.show_aln(Alignment_v2.calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i+1],delta))
        check_aln1 += str(aln.seq1[start_seq1:delta])
        check_aln2 += str(aln.seq2[start_seq2:i+1])
        print ""
      
      elif i == tp[-1]:
        print count * delta, "...", count * delta + delta - 1
        print tp[count - 1] + 1, "...", i
        aln = Alignment_v2.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        # print Alignment_v2.show_aln(Alignment_v2.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta))
        print Alignment_v2.show_aln(aln.seq1, aln.seq2)
        check_aln1 += str(aln.seq1[count*delta:count*delta+delta])
        check_aln2 += str(aln.seq2[tp[count-1]+1:i+1])
        print ""
        count += 1
        print count * delta, "...", len(seq1) - 1
        print tp[count - 1] + 1, "...", len(seq2) - 1
        aln = Alignment_v2.calculate_alignment(seq1[count*delta:len(seq1)],seq2[tp[count-1]+1:len(seq2)],delta)
        # print Alignment_v2.show_aln(Alignment_v2.calculate_alignment(seq1[count*delta:len(seq1)],seq2[tp[count-1]+1:len(seq2)],delta))
        print Alignment_v2.show_aln(aln.seq1, aln.seq2)
        check_aln1 += str(aln.seq1[count*delta:len(seq1)])
        check_aln2 += str(aln.seq2[tp[count-1]+1:len(seq2)])
        break
      else:
        print count * delta, "...", count * delta + delta - 1
        print tp[count - 1] + 1, "...", i
        aln = Alignment_v2.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        # print Alignment_v2.show_aln(Alignment_v2.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta))
        print Alignment_v2.show_aln(aln.seq1, aln.seq2)
        check_aln1 += str(aln.seq1[count*delta:count*delta+delta])
        check_aln2 += str(aln.seq2[tp[count-1]+1:i+1])
        print ""
        count += 1
    
    return True
  
  # store TracePointAlignment to file
    def store_aln(self, output, mode):
  
      with open('alignment_compressed.txt', mode) as file_:
  
        # ';' for easy splitting
        for item in output:
          file_.write("%s;" % item)

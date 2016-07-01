# -*- coding: utf-8 -*-

import Alignment

import re
import math
import sys

class TracePointAlignment(object):

  def __init__(self, seq1, seq2, delta, start_seq1, start_seq2):
    self.seq1 = seq1
    self.seq2 = seq2
    self.delta = delta
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2
    # TODO geht das so?
    self.tp = None

  # extract TracePoints from CIGAR-String
  def encode_cigar(self, cigar):

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

    print "Trace Points:", tp
    self.tp = tp

    # kein TPAln...
    # tp_aln = TracePoint_v3.TracePointAlignment(self.seq1, self.seq2, self.delta, tp=tp, start_seq1=self.start_seq1, start_seq2=self.start_seq2)
    # return tp_aln

  #number of chars and '-'s in sequence  
  def count_indels_letters(self, seq):

    letter_count = 0
    indel_count = 0
    for i in seq:
      if i.isalpha():
        letter_count += 1
      elif i == '-':
        indel_count += 1

    return [letter_count, indel_count]

  # extract TracePoints from alignment
  def encode(self, alignment):
  
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

    print "# Trace Points:", tp, "\n"
  
    return TP_alignment

  # create new intervals from TracePoints and calculate new alignment
  def decode(self, seq1, seq2, delta, tp, start_seq1, start_seq2):
    count = 0
    new_seq1 = new_seq2 = ""
    for i in tp:
      if i == tp[0]:
        count += 1
        aln = Alignment.calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
      elif i == tp[-1]:
        aln = Alignment.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        count += 1
        aln = Alignment.calculate_alignment(seq1[count*delta:len(seq1)],seq2[tp[count-1]+1:len(seq2)],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        break
      else:
        aln = Alignment.calculate_alignment(seq1[count*delta:count*delta+delta],seq2[tp[count-1]+1:i+1],delta)
        new_seq1 += str(aln.seq1)
        new_seq2 += str(aln.seq2)
        count += 1

    
    print "# Konkateniertes Alignment:"
    print Alignment_v3.show_aln(new_seq1, new_seq2)
    
    return True
  
  # TODO implementation not finished yet
  # store TracePointAlignment to file
    def store_aln(self, output, mode):
  
      with open('alignment_compressed.txt', mode) as file_:
  
        # ';' for easy splitting
        for item in output:
          file_.write("%s;" % item)

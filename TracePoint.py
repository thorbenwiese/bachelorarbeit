import Alignment

import re
import math
import sys

class TracePointAlignment(object):

  def __init__(self, seq1, seq2, delta, start_seq1, end_seq1, start_seq2, end_seq2):
    self.seq1 = seq1
    self.seq2 = seq2
    self.delta = delta
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2
    self.end_seq1 = end_seq1
    self.end_seq2 = end_seq2
    # is set by encode function
    self.tp = None

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

  # extract TracePoints from CIGAR-String
  def encode(self, cigar):

    count = count1 = count2 = cig_count = interval_count = 0

    # dynamic calculation of interval size
    dynamic = max(1,int(math.ceil(self.start_seq1/self.delta)))

    # adjustment of interval length
    interval_count = min(int(math.ceil(float(len(self.seq1)) / self.delta)), 
                         int(math.ceil(float(len(self.seq2)) / self.delta)))

    # intervals
    # initialized with 0s
    intervals = [0] * interval_count

    for i in range(0, interval_count):
      # first interval
      if i == 0:
        intervals[i] = self.start_seq1, dynamic * self.delta - 1
      # last interval
      elif i == interval_count - 1:
        intervals[i] = (dynamic + interval_count - 2) * self.delta, self.end_seq1-1
      # other intervals
      else:
        intervals[i] = (dynamic + i - 1) * self.delta, (dynamic + i) * self.delta - 1

    # create pattern for CIGAR-String
    cigar_pattern = re.compile(r"\d+[MIDNSHP=j]{1}")

    # search cigar for pattern
    tp = []
    for j in cigar_pattern.findall(cigar):
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
          tp.append(count2 - 1 + self.start_seq2)
          if count != len(intervals)-1:
            count += 1

    self.tp = tp

  # TODO macht gar nichts :D
  # create new intervals from TracePoints and calculate new alignment
  def decode(self, seq1, seq2, delta, tp, start_seq1, start_seq2):

    new_seq1 = new_seq2 = ""
    
    aln = Alignment.Alignment(seq1, seq2, start_seq1, start_seq2)

    for i in range(0,len(tp)):
      
      if i == 0:

        new_seq1 += str(seq1[0:delta])
        new_seq2 += str(seq2[0:tp[i] + 1])
      
      elif i == len(tp) - 1:

        new_seq1 += str(seq1[i*delta:(i+1)*delta])
        new_seq2 += str(seq2[tp[i-1]+1:tp[i]+1])
 
        new_seq1 += str(seq1[(i+1)*delta:len(seq1)])
        new_seq2 += str(seq2[tp[i]+1:len(seq2)])

      else:
        
        new_seq1 += str(seq1[i*delta:(i+1)*delta])
        new_seq2 += str(seq2[tp[i-1]+1:tp[i] + 1])

  # store TracePointAlignment to file
  def store_tp_aln(self):
  
    with open('aln_file.txt', 'a') as file_:
      
      # mit join?
      #file_.write(";".join((self.delta, self.start_seq1, self.end_seq1,
      #                     self.start_seq2, self.end_seq2, self.tp)))
      file_.write("%d;%d;%d;%d;%d;%s\n" % (self.delta, self.start_seq1, 
                  self.end_seq1, self.start_seq2, self.end_seq2, self.tp))

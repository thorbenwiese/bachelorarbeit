import Alignment
import Cigar

import re
import math
import sys

# pattern for CIGAR-String
cigar_pattern = re.compile(r"\d+[MID]{1}")

class TracePointAlignment(object):

  def __init__(self, seq1, seq2, start_seq1, end_seq1, start_seq2, end_seq2, 
               delta, cigar = None):
    self.seq1 = seq1
    self.seq2 = seq2
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2
    self.end_seq1 = end_seq1
    self.end_seq2 = end_seq2
    self.delta = delta
    self.cigar = cigar
    if cigar != None:
      self.tp = self.encode()

  # extract TracePoints from CIGAR-String
  def encode(self):

    count = count1 = count2 = cig_count = interval_count = 0

    # calculation of interval size
    itv_size = max(1,int(math.ceil(self.start_seq1/self.delta)))

    # adjustment of interval length
    interval_count = min(int(math.ceil(float(len(self.seq1)) / self.delta)), 
                         int(math.ceil(float(len(self.seq2)) / self.delta)))

    # intervals
    intervals = [0] * interval_count

    for i in range(0, interval_count):
      
      if i == 0: 
        begin = self.start_seq1

      else:
        begin = (itv_size+ i - 1) * self.delta

      if i == interval_count - 1:
        end = self.end_seq1 - 1

      else:
        end = (itv_size+i) * self.delta - 1
      
      intervals[i] = begin, end

    # search cigar for pattern
    tp = []
    for j in cigar_pattern.findall(self.cigar):

      cig_count = int(j[:-1])
      cig_symbol = j[1]
      assert cig_symbol not in ['N','S','H','P'], \
        "CIGAR-Symbol is not in ['M','I','D']"
   
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
          if count != len(intervals) - 1:
            count += 1
    assert tp, "TracePoint Array from encode function is empty."

    return tp

  # create new intervals from TracePoints and calculate new alignment
  def decode(self, tp):

    assert self.seq1, "First sequence for decode function is empty."
    assert self.seq2, "Second sequence for decode function is empty."
    assert self.delta > 0, "Delta for decode function is <= 0."
    assert tp, "TracePoint Array for decode function is empty."
    assert self.start_seq1 >= 0, \
      "Starting position for first sequence in decode function is < 0."
    assert self.start_seq2 >= 0, \
      "Starting position for second sequence in decode function is < 0."
    assert self.end_seq1 > 0, \
      "End position for first sequence in decode function is <= 0."
    assert self.end_seq2 > 0, \
      "End position for second sequence in decode function is <= 0."

    # calculate CIGAR of intervals
    cigar = ""
    
    aln = Alignment.Alignment(self.seq1, self.seq2, self.start_seq1, 
                              self.end_seq1,self.start_seq2,self.end_seq2)

    for i in range(0,len(tp)):
      
      if i == 0:

        cigar = aln.calc_cigar(self.seq1[0:self.delta],self.seq2[0:tp[i]+1])
      
      elif i == len(tp) - 1:
 
        cigar += aln.calc_cigar(self.seq1[i*self.delta:len(self.seq1)],
                                self.seq2[tp[i-1]+1:len(self.seq2)])

      else:
        
        cigar += aln.calc_cigar(self.seq1[i*self.delta:(i+1)*self.delta],
                                      self.seq2[tp[i-1]+1:tp[i] + 1])

    print "CIGAR:", cigar
    cigar = Cigar.combine_cigar(cigar)
    print "CIGAR2:", cigar

    return cigar

  # store TracePointAlignment to file
  def store_tp_aln(self,mode):
  
    with open('aln_file.txt', mode) as file_:
      
      file_.write("%d;%d;%d;%d;%d;%s\n" % (self.delta, self.start_seq1, 
                  self.end_seq1, self.start_seq2, self.end_seq2, self.tp))

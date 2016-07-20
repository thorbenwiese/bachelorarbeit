import Alignment

import re
import math
import sys

class TracePointAlignment(object):

  def __init__(self, seq1, seq2, start_seq1, end_seq1, start_seq2, end_seq2, 
               delta, cigar = None):
    self.seq1 = seq1
    self.seq2 = seq2
    self.delta = delta
    self.start_seq1 = start_seq1
    self.start_seq2 = start_seq2
    self.end_seq1 = end_seq1
    self.end_seq2 = end_seq2
    self.cigar = cigar
    if cigar != None:
      self.tp = self.encode()

  # extract TracePoints from CIGAR-String
  def encode(self):

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
        intervals[i] = (dynamic+interval_count-2) * self.delta, self.end_seq1-1
      # other intervals
      else:
        intervals[i] = (dynamic + i-1) * self.delta, (dynamic+i) * self.delta-1

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
          tp.append(count2 - 1 + self.start_seq2)
          if count != len(intervals)-1:
            count += 1

    assert tp, "TracePoint Array from encode function is empty."

    return tp

  # create new intervals from TracePoints and calculate new alignment
  def decode(self, tp):

    assert self.seq1, "First sequence for decode function is empty."
    assert self.seq2, "Second sequence for decode function is empty."
    assert self.delta > 0, "Delta for decode function is <= 0."
    assert tp, "TracePoint Array for decode function is empty."
    assert self.start_seq1 >= 0, "Starting position for first sequence in decode function is < 0."
    assert self.start_seq2 >= 0, "Starting position for second sequence in decode function is < 0."
    assert self.end_seq1 > 0, "End position for first sequence in decode function is <= 0."
    assert self.end_seq2 > 0, "End position for second sequence in decode function is <= 0."

    # calculate CIGAR of intervals
    cigar = ""
    
    aln = Alignment.Alignment(self.seq1, self.seq2, self.start_seq1, 
                              self.end_seq1,self.start_seq2,self.end_seq2)

    for i in range(0,len(tp)):
      
      if i == 0:

        aln1, aln2 = aln.calculate(self.seq1[0:self.delta],self.seq2[0:tp[i]+1])
        cigar += aln.calc_cigar(aln1, aln2)
      
      elif i == len(tp) - 1:
 
        aln1, aln2 = aln.calculate(self.seq1[i*self.delta:len(self.seq1)],
                                   self.seq2[tp[i-1]+1:len(self.seq2)])
        cigar += aln.calc_cigar(aln1, aln2)

      else:
        
        aln1, aln2 = aln.calculate(self.seq1[i*self.delta:(i+1)*self.delta],
                                   self.seq2[tp[i-1]+1:tp[i] + 1])
        cigar += aln.calc_cigar(aln1, aln2)

    # calculate aln_seq with CIGAR

    cig_count = tmp1 = tmp2 = count = 0
    aln1 = aln2 = ""

    #new pattern for CIGAR-String, format: Amount + 1 char from {M,I,D,N,S,H,P}
    cigar_pattern = re.compile(r"\d+[MIDNSHP=X]{1}")

    #in cigar nach pattern suchen
    for element in cigar_pattern.findall(cigar):
      count+=1

      tmp1 += cig_count
      tmp2 += cig_count
      cig_count = int(element[:-1])
      cig_symbol = element[-1]

      if cig_symbol == 'M':
        aln1 += str(self.seq1[tmp1:tmp1+cig_count])
        aln2 += str(self.seq2[tmp2:tmp2+cig_count])

      elif cig_symbol == 'I':
        aln1 += str(self.seq1[tmp1:tmp1+cig_count])
        aln2 += str("-" * cig_count)
        tmp2-=cig_count

      elif cig_symbol == 'D':
        aln1 += str("-" * cig_count)
        aln2 += str(self.seq2[tmp2:tmp2+cig_count])
        tmp1 -= cig_count

    aln.show_aln(aln1, aln2)


  # store TracePointAlignment to file
  def store_tp_aln(self,mode):
  
    with open('aln_file.txt', mode) as file_:
      
      file_.write("%d;%d;%d;%d;%d;%s\n" % (self.delta, self.start_seq1, 
                  self.end_seq1, self.start_seq2, self.end_seq2, self.tp))

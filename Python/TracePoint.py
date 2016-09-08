import Alignment
import Cigar_Pattern
import math

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

    # calculation of interval size
    # p is a factor to dynamically adjust the interval borders
    # if the sequence starts at pos 0 then p should be 1
    p = max(1,int(math.ceil(self.start_seq1/self.delta)))

    # number of intervals
    tau = int(math.ceil(float(self.end_seq1) / self.delta) - math.floor(
                              self.start_seq1 / self.delta))

    # List of last indices of intervals in u initialized with 0
    # tau - 1 because the last interval has no Trace Point
    u_tp = [0] * (tau - 1)

    # tau - 1 because last interval has no Trace Point
    for q in range(0, len(u_tp)):
      u_tp[q] = (p + q) * self.delta - 1

    v_tp = []
    num_chars_in_u = num_chars_in_v = count = 0

    for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(self.cigar):

      assert cig_symbol in ['M', 'I', 'D'], \
        "CIGAR-Symbol is not in ['M','I','D']"
   
      for i in range(0,cig_count):
        if cig_symbol == 'I':
          num_chars_in_v += 1
        elif cig_symbol == 'D':
          num_chars_in_u += 1
        else:
          num_chars_in_u += 1
          num_chars_in_v += 1

        # count until the end but ignore end of last interval as TracePoint
        if num_chars_in_u == u_tp[count]:
          v_tp.append(num_chars_in_v)

          # do not increment count if the last element in u_tp is reached
          if count == len(u_tp) - 1:
            assert v_tp, "TracePoint Array from encode function is empty."
            return v_tp

          else:
            count += 1



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

        cigar = aln.calc_cigar(self.seq1[0:self.delta], self.seq2[0:tp[i] + 1])
      
      elif i == len(tp) - 1:
 
        cigar += aln.calc_cigar(self.seq1[i*self.delta:len(self.seq1)],
                                self.seq2[tp[i - 1] + 1:len(self.seq2)])

      else:
        
        cigar += aln.calc_cigar(self.seq1[i * self.delta:(i + 1) * self.delta],
                                self.seq2[tp[i - 1] + 1:tp[i] + 1])

    cigar = Cigar_Pattern.combine_cigar(cigar)

    return cigar

  # store TracePointAlignment to file
  def store_tp_aln(self,mode):
  
    with open('aln_file.txt', mode) as file_:
      
      file_.write("%d;%d;%d;%d;%d;%s\n" % (self.delta, self.start_seq1, 
                  self.end_seq1, self.start_seq2, self.end_seq2, self.tp))

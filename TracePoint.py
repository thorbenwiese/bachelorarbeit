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
    # is set by encode function
    self.tp = None

  # TODO in one loop without raw_cigar -> less memory
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



    """  
      # deletion 
      if aln_seq1[i] == '-':
        while aln_seq1[i] == '-':
          if i == seqlen:
            break;
          else:
            i += 1
        count = i - count
        cigar += "%d%s" % (count+1, 'D')
        print "COUNT", count, "I", i, "CIGAR", cigar
      # insertion
      elif aln_seq2[i] == '-':
        while aln_seq2[i] == '-':
          if i == seqlen:
            break;
          else:
            i += 1
        count = i - count
        cigar += "%d%s" % (count+1, 'I')
        print "COUNT", count, "I", i, "CIGAR", cigar
      # match
      elif aln_seq1[i] == aln_seq2[i] or aln_seq1[i] != aln_seq2[i]:
        while aln_seq1[i] == aln_seq2[i] or aln_seq1[i] != aln_seq2[i]:
          if i == seqlen:
            break;
          else:
            i += 1
          print "IIIII:    ",i
        count = i - count
        cigar += "%d%s" % (count+1, 'M')
        print "COUNT", count, "I", i, "CIGAR", cigar
    print "CIGAR:", cigar
    ###  
      if aln_seq1[i+1] == aln_seq2[i+1]:
          count += 1
      elif aln_seq1[i] == '-':
        raw_cigar += 'D'
      elif aln_seq2[i] == '-':
        raw_cigar += 'I'

    ####

    raw_cigar = ""
    for i in range(0, seqlen):
      if aln_seq1[i] == aln_seq2[i]:
        raw_cigar += 'M'
      elif aln_seq1[i] == '-':
        raw_cigar += 'D'
      elif aln_seq2[i] == '-':
        raw_cigar += 'I'

    cigar = ""

    count = 1
    for j in range(0, len(raw_cigar)):

      if j == len(raw_cigar)-1:
        cigar += "%d%s" % (count, raw_cigar[j-1])
        break
      elif raw_cigar[j] == raw_cigar[j+1]:
        count += 1
      else:
        cigar += "%d%s" % (count, raw_cigar[j])
        count = 1
    """
    return cigar

  # extract TracePoints from CIGAR-String
  def encode_cigar(self, cigar):

    count = count1 = count2 = cig_count = interval_count = 0

    end_seq1 = self.start_seq1 + len(self.seq1) - 1
    end_seq2 = self.start_seq2 + len(self.seq2) - 1

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
          count += 1

    self.tp = tp

  #number of chars and '-'s in sequence  
  def count_indels_letters(self, seq):

    letter_count = 0
    indel_count = 0

    for i in range(0,len(seq)):
      if seq[i].isalpha():
        letter_count += 1
      elif seq[i] == '-':
        indel_count += 1

    return [letter_count, indel_count]

  # extract TracePoints from alignment
  def encode(self):
  
    indels_in_seq2 = all_chars_in_seq1 = interval_count = 0

    # number of chars in sequences
    letters_in_seq1 = self.count_indels_letters(self.seq1)[0]
    letters_in_seq2 = self.count_indels_letters(self.seq2)[0]

    # Trace Point Array
    tp = []

    end_seq1 = letters_in_seq1 - 1
    end_seq2 = letters_in_seq2 - 1

    # dynamic calculation of interval size
    dynamic = max(1,int(math.ceil(self.start_seq1/self.delta)))

    # adjustment of interval length
    interval_count = min(int(math.ceil(float(len(self.seq1)) / self.delta)), int(math.ceil(float(len(self.seq2)) / self.delta)))

    # intervals
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

    start  = tmp = 0

    for i in range(0,len(intervals)):
    
      # ignore last interval
      if i == len(intervals) - 2:
        break
      else:

        # Anzahl der Buchstaben an den Enden der Sequenzen
        rest1 = self.count_indels_letters(self.seq1[all_chars_in_seq1:end_seq1])[0]
        rest2 = self.count_indels_letters(self.seq2[all_chars_in_seq1:end_seq2])[0]

        # break condition
        if rest1 != 0 and rest2 != 0:
          # first delta letters in alignment.seq1
          while self.count_indels_letters(self.seq1[start:all_chars_in_seq1])[0] != self.delta:
            all_chars_in_seq1 += 1

        # count gaps in seq2
        indel_count = self.count_indels_letters(self.seq2[start:all_chars_in_seq1])[1]
        indels_in_seq2 += indel_count
        tp.append(all_chars_in_seq1 - indels_in_seq2 - 1 + self.start_seq2)
        tmp = start
        start = all_chars_in_seq1
        all_chars_in_seq1 = tmp

    self.tp = tp

    self.store_tp_aln()

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
      
      file_.write("%d;%d;%d;%s\n" % (self.delta, self.start_seq1, self.start_seq2, self.tp))


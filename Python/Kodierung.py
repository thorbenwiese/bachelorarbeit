#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tp_calc
import Alignment
import TracePoint
import Cigar_Pattern
import math
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
import matplotlib.pyplot as plt
import collections
import argparse
import time
import sys
reload(sys)
sys.setdefaultencoding('utf8')

from heapq import heappush, heappop, heapify
from collections import defaultdict

def bucket(bit_sum, base):

  bucket = []
  for element in bit_sum:
    r = int(base * round(float(element)/base))
    if r == 0:
      bucket.append(element)
    else:
      bucket.append(r)

  return bucket

def encode(symb2freq):
    """Huffman encode the given dict mapping symbols to weights"""
    heap = [[wt, [sym, ""]] for sym, wt in symb2freq.items()]
    heapify(heap)
    while len(heap) > 1:
        lo = heappop(heap)
        hi = heappop(heap)
        for pair in lo[1:]:
            pair[1] = '0' + pair[1]
        for pair in hi[1:]:
            pair[1] = '1' + pair[1]
        heappush(heap, [lo[0] + hi[0]] + lo[1:] + hi[1:])
    return sorted(heappop(heap)[1:], key=lambda p: (len(p[-1]), p))
 
def huffman(code):
  symb2freq = defaultdict(int)
  for ch in code:
    symb2freq[ch] += 1
  huff = encode(symb2freq)
  bits = 0
  canon = []
  # print "Symbol\tWeight\tHuffman Code"
  for p in huff:
    bits += len(p[1])
    it = (p[0],p[1],len(p[1]))
    canon.append(it)
    # print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
  bits += canonical(canon)
  return bits 

def canonical(huf):
  #for a in huf:print a
  minsize=huf[0][2]
  co = '0'
  if len(co)<minsize:co='0'*(minsize-len(co))+co
  ne = [(huf[0][0],co,len(co))]
  code = 0
  for i in range(1,len(huf)):
    code = (code+1) << huf[i][2] - huf[i-1][2]
    co=bin(code)[2:]
    l = len(co)
    if l<minsize:co='0'*(minsize-l)+co
    ne.append((huf[i][0],co,len(co)))
  #for a in ne:
  #  print a
  ls1 = []
  ls2 = []
  ls3 = []
  for a in ne:
    ls1.append(a[2]) # len(bits)
    ls2.append(a[1]) # bits
    ls3.append(a[0]) # symbol

  c = collections.Counter(ls1)
  vals = c.values()

  for i in range(0,len(ls2)):
    if len(ls2[i]) != i+1:
      vals.insert(i,0)
  leng = len(ls2)-1
  for i in range(leng,0,-1):
    if vals[i] == 0:
      del vals[i]
    elif vals[i] != 0:
      break
  size_header = 0
  for element in ls1:
    size_header+=(element+1)
  size_header+=(len(ls3))
  
  return size_header

def unary(code):
  symb2freq = defaultdict(int)
  for ch in code:
    symb2freq[ch] += 1
  huff = encode(symb2freq)
  bits = 0
  i = 1
  order = []
  for p in huff:
    order.append(symb2freq[p[0]])
  order.sort(reverse=True) 
  for i in range(0,len(order)):
    bits += (i+1)*order[i]
  return bits 

def entropy(amount, random_length, error_rate, alphabet, delta):
  start_seq1 = start_seq2 = 0
  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  diff_ent_sum = []
  cig_ent_sum = []
  for i in range(0, len(random_seq_list), 2):
    end_seq1 = len(random_seq_list[i])
    end_seq2 = len(random_seq_list[i+1]) 
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                    start_seq1, end_seq1, start_seq2, end_seq2)

    cigar = aln.calc_cigar(aln.seq1, aln.seq2)
    tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1, 
                                            end_seq1, start_seq2, end_seq2, 
                                            delta, cigar)
    mcount = dcount = icount = 0
    cig_count_sum = []
    cig_sym_sum = []
    for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
      cig_sym_sum.append(cig_symbol)
      cig_count_sum.append(cig_count)
    counter1 = collections.Counter(cig_count_sum)    
    counter2 = collections.Counter(cig_sym_sum)    

    ent1 = ent2 = 0
    for value in counter1.values():
      px = float(value)/sum(counter1.values())
      ent1 += px * math.log(px,2)

    for value in counter2.values():
      px = float(value)/sum(counter2.values())
      ent2 += px * math.log(px,2)
    cig_ent_sum.append(-ent1 * len(cig_count_sum) + (-ent2 * len(cig_sym_sum)))

    # store differences between Trace Points and Delta Value in List
    TP = [delta, tp_aln.tp[0]]
    for j in range(1,len(tp_aln.tp)):
      TP.append(tp_aln.tp[j]-tp_aln.tp[j-1])

    counter = collections.Counter(TP)
    #counter = collections.Counter(bucket(TP))
    ent = 0
    for value in counter.values():
      px = float(value)/sum(counter.values())
      ent += px * math.log(px,2)
    diff_ent_sum.append(-ent*len(TP))

  cig_counter = collections.Counter(bucket(cig_ent_sum,10))
  diff_counter = collections.Counter(bucket(diff_ent_sum,3))

  cig_mean = diff_mean = 0

  plt.figure(1)
  plt.ylabel('Anzahl Sequenzpaare')
  plt.xlabel('Entropie in Bit')
  plt.axis([0, max(cig_counter.keys())*1.2, 0, 
            cig_counter.most_common(1)[0][1] * 1.2])
  plt.plot(cig_counter.keys(), cig_counter.values(), 
           'bo', label="CIGAR Entropy")
  cig_mean = float(sum(cig_ent_sum))/len(cig_ent_sum)
  y = list(range(cig_counter.most_common(1)[0][1]))
  x = [cig_mean] * cig_counter.most_common(1)[0][1]

  print "Cig Entropy Mean:", cig_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/cig_entropy")

  plt.figure(2)
  plt.ylabel('Anzahl Sequenzpaare')
  plt.xlabel('Entropie in Bit')
  plt.axis([0, max(diff_counter.keys())*1.2, 0, 
            diff_counter.most_common(1)[0][1] * 1.2])
  plt.plot(diff_counter.keys(), diff_counter.values(), 
           'bo', label="Differences Entropy")
  diff_mean = float(sum(diff_ent_sum))/len(diff_ent_sum)
  y = list(range(diff_counter.most_common(1)[0][1]))
  x = [diff_mean] * diff_counter.most_common(1)[0][1]
  print "Diff Entropy Mean:", diff_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/diff_entropy")

  print "cig/diff:", float(cig_mean)/diff_mean
  print "diff/cig:", float(diff_mean)/cig_mean
  print cig_counter
  print diff_counter
  print [int(i) for i in cig_counter.keys()]
  print diff_counter.keys()

  plt.show()



def cigar_kodierung(ciglist, mode):
  bit_sum = []

  for cigar in ciglist:
    if mode == "huffman":
      cig_count_sum = []
      bits = 2
      for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
        bits += 1
        #if cig_symbol == 'M':
          # bit_sum.append(1)
         # bits += 1
        #else:
          # bit_sum.append(2)
          #bits += 2
        cig_count_sum.append(cig_count)
      bit_sum.append(huffman(cig_count_sum) + bits)
    elif mode == "binary":
      # naiive binary coding for cigar
      bin_count = 0

      for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
        # cig_count_sum.append(cig_count)
        bin_count += 2
      # bit_sum.append((int(math.ceil(math.log(len(cig_count_sum),2))) * len(
      #                cig_count_sum)) + 2 * len(cig_count_sum))
      bit_sum.append(int(math.ceil(math.log(bin_count,2))) * bin_count)
      
    else:
      # unary coding
      m_count = i_count = d_count = add = 0
      cig_count_sum = []
      for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
        if cig_symbol == 'M':
          m_count += 1
        elif cig_symbol == 'D':
          d_count += 1
        else:
          i_count += 1
        cig_count_sum.append(cig_count)
        if d_count < i_count:
          add = m_count + i_count * 2 + d_count * 3
        else:
          add = m_count + d_count * 2 + i_count * 3
      bit_sum.append(unary(cig_count_sum) + add)

  return bit_sum

def tp_kodierung(TP, mode):

  bit_sum = []
  if mode == "binary":
    # naiive binary coding for differences
    bit_sum.append(int(math.ceil(math.log(len(TP),2)))*len(TP))
  elif mode == "unary":
    # unary coding for differences
    bit_sum.append(unary(TP))
  else:
    # huffman coding for differences
    bit_sum.append(huffman(TP))
  return bit_sum

def multiplot(bs1, bs2, bs3, buck1, buck2, buck3, method, t):

  
  counter1 = collections.Counter(bucket(bs1, buck1))
  counter2 = collections.Counter(bucket(bs2, buck2))
  counter3 = collections.Counter(bucket(bs3, buck3))

  bin_mean = una_mean = huf_mean = 0

  if method == "tracepoint":

    plt.figure(1)
    plt.ylabel('Anzahl der Sequenzpaare')
    plt.xlabel('Größe der Kodierung in Bit')

    maxkey = max(max(counter2.keys()),max(counter3.keys()))
    maxend = max(counter2.most_common(1)[0][1],
                 counter3.most_common(1)[0][1])
    plt.axis([0, maxkey * 1.2, 0, maxend * 1.2])
  
    plt.plot(counter2.keys(), counter2.values(), 'co', 
             label="Unär")
    plt.plot(counter3.keys(), counter3.values(), 'ro', 
             label="Huffman")
    plt.legend(loc='upper right')

  else: # cigar

    plt.figure(1)
    plt.ylabel('Anzahl der Sequenzpaare')
    plt.xlabel('Größe der Kodierung in Bit')

    maxkey = max(max(counter1.keys()),max(counter2.keys()),max(counter3.keys()))
    maxend = max(counter1.most_common(1)[0][1],counter2.most_common(1)[0][1],
                 counter3.most_common(1)[0][1])
    #maxkey = max(max(counter1.keys()),max(counter2.keys()))
    #maxend = max(counter1.most_common(1)[0][1],counter2.most_common(1)[0][1])
    plt.axis([0, maxkey * 1.2, 0, maxend * 1.2])
  
    plt.plot(counter1.keys(), counter1.values(),'bo', 
             label="Binär")
    plt.plot(counter2.keys(), counter2.values(), 'co', 
             label="Unär")
    plt.plot(counter3.keys(), counter3.values(), 'ro', 
             label="Huffman")
    plt.legend(loc='upper right')

    """
    plt.figure(2)
    plt.ylabel('Anzahl der Sequenzpaare')
    plt.xlabel('Größe der Kodierung in Bit')
    plt.axis([0, max(counter3.keys()) * 1.2, 
              0, counter3.most_common(1)[0][1] * 1.2])

    plt.plot(counter3.keys(), counter3.values(), 'ro', 
           label="Huffman")
    plt.legend(loc='upper right')
    """

  bin_mean = float(sum(bs1))/len(bs1)
  una_mean = float(sum(bs2))/len(bs2)
  huf_mean = float(sum(bs3))/len(bs3)

  print ""
  print "Binärer Counter:\t", counter1
  print ""
  print "Unärer Counter:\t\t", counter2
  print "Unäre BITSUM:\t\t", bs2
  print ""
  print "Huffman Counter:\t", counter3
  print "##########################"
  print "Binäre Keys:\t", counter1.keys()
  print ""
  print "Unäre Keys:\t", counter2.keys()
  print ""
  print "Huffman Keys:\t", counter3.keys()
  print ""
  print "Bin/Una:",float(bin_mean)/una_mean
  print "Una/Bin:",float(una_mean)/bin_mean
  print ""

  print "Bin/Huf:",float(bin_mean)/huf_mean
  print "Huf/Bin:",float(huf_mean)/bin_mean
  print ""

  print "Huf/Una:",float(huf_mean)/una_mean
  print "Una/Huf:",float(una_mean)/huf_mean
  print ""


  print "MAX Bin:", max(counter1.keys())
  print "MIN Bin:", min(counter1.keys())
  print "Binary Mean:", bin_mean
  print ""

  print "MAX Una:", max(counter2.keys())
  print "MIN Una:", min(counter2.keys())
  print "Unary Mean:", una_mean
  print ""

  print "MAX Huf:", max(counter3.keys())
  print "MIN Huf:", min(counter3.keys())
  print "Huffman Mean:", huf_mean
  print ""

  print "Calculation complete.\nClock time: %.2f seconds." % (time.clock() - t)
  plt.show()






def main():

 #  METHOD = "cigar" 
  METHOD = "tracepoint"
  # METHOD = "entropy"
  start_seq1 = start_seq2 = 0

  t = time.clock()
  print "Läuft..."
  amount = 100
  length = 5000
  delta = 100
  err_rate = 0.15
  bs1 = []
  bs2 = []
  bs3 = []
  TP = []

  random_seq_list = tp_calc.random_sequences(amount, length, err_rate, "acgt")

  print METHOD
  
  ciglist = []
  for i in range(0, len(random_seq_list), 2):
    print i/2 + 1, "von", len(random_seq_list)/2
    end_seq1 = len(random_seq_list[i])
    end_seq2 = len(random_seq_list[i + 1])
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                              start_seq1, end_seq1, start_seq2, end_seq2)

    if METHOD == "cigar":

      # CIGAR

      ciglist.append(aln.calc_cigar(aln.seq1, aln.seq2))

      if i == len(random_seq_list) - 2:
        
        print "LÄNGE:", len(ciglist)

        bs1 = cigar_kodierung(ciglist, "binary")
        bs2 = cigar_kodierung(ciglist, "unary")
        bs3 = cigar_kodierung(ciglist, "huffman")

        multiplot(bs1, bs2, bs3, 140, 100, 25, METHOD, t)
  
    elif METHOD == "tracepoint":

      # Differenzen

      cigar = aln.calc_cigar(aln.seq1, aln.seq2)
      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1,
                                end_seq1, start_seq2, end_seq2, delta, cigar)

      if i == len(random_seq_list) - 2:
        multiplot(bs1, bs2, bs3, 1, 500, 25, METHOD, t)
      else:
        TP.append(delta)
        TP.append(tp_aln.tp[0])
        for j in range(1, len(tp_aln.tp)):
          TP.append(tp_aln.tp[j] - tp_aln.tp[j - 1])

        for a in tp_kodierung(TP, "binary"): bs1.append(a)
        for b in tp_kodierung(TP, "unary"): bs2.append(b)
        for c in tp_kodierung(TP, "tracepoint"): bs3.append(c)
  
    else:

      # Entropy
      entropy(amount,length,err_rate,"acgt",delta)

if __name__ == "__main__":
  main()

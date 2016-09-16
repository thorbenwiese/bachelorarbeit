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

def create_cigars_or_tp(amount, random_length, error_rate, alphabet, delta, method):
  start_seq1 = start_seq2 = 0
  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  if method == "cigar":
    ciglist = []
    for i in range(0, len(random_seq_list), 2):
      end_seq1 = len(random_seq_list[i])
      end_seq2 = len(random_seq_list[i+1]) 
      aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                      start_seq1, end_seq1, start_seq2, end_seq2)

      ciglist.append(aln.calc_cigar(aln.seq1, aln.seq2))

    return [ciglist, delta]
    #return ["4M1I1M1I1M1I1M2D1M1D1M1I8M1D7M1D5M1I4M"]

  else:
    
    for i in range(0, len(random_seq_list), 2):
      end_seq1 = len(random_seq_list[i])
      end_seq2 = len(random_seq_list[i+1]) 
      aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                      start_seq1, end_seq1, start_seq2, end_seq2)
      cigar = aln.calc_cigar(aln.seq1, aln.seq2)
      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, 
                      start_seq1, end_seq1, start_seq2, end_seq2, delta, cigar)

      # store differences between Trace Points and Delta Value in List
      TP = [delta, tp_aln.tp[0]]
      for j in range(1,len(tp_aln.tp)):
        TP.append(tp_aln.tp[j]-tp_aln.tp[j-1])
    return[TP,delta]

def kodierung(method, mode, delta, ciglist=None, TP=None):

  bit_sum = []
  cig_count_sum = []

  for cigar in ciglist:
    # huffman coding for cigar
    if method == "cigar":
      if mode == "huffman":
        for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
          if cig_symbol == 'M':
            bit_sum.append(1)
          else:
            bit_sum.append(2)
          cig_count_sum.append(cig_count)
        bit_sum.append(huffman(cig_count_sum))
        cig_count_sum = []
      elif mode == "binary":
        # naiive binary coding for cigar
        for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
          cig_count_sum.append(cig_count)
        bit_sum.append(int(math.ceil(math.log(len(cig_count_sum),2)))*len(
                           cig_count_sum))
        bit_sum.append(2 * len(cig_count_sum))
        cig_count_sum = []
      else:
        # unary coding
        m_count = i_count = d_count = 0
        for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
          if cig_symbol == 'M':
            m_count += 1
          elif cig_symbol == 'D':
            d_count += 1
          else:
            i_count += 1
          cig_count_sum.append(cig_count)
        bit_sum.append(unary(cig_count_sum))
        if d_count < i_count:
          bit_sum.append(m_count + i_count * 2 + d_count * 3)
        else:
          bit_sum.append(m_count + d_count * 2 + i_count * 3)
        cig_count_sum = []
    else:
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

def bucket(bit_sum, base):

  bucket = []
  for element in bit_sum:
    r = int(base * round(float(element)/base))
    if r == 0:
      bucket.append(element)
    else:
      bucket.append(r)

  return bucket

def multiplot(bs1, bs2, bs3, t, buck, method):

  counter1 = collections.Counter(bucket(bs1,buck))
  counter2 = collections.Counter(bucket(bs2,buck))
  counter3 = collections.Counter(bucket(bs3,buck))

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

  else:
    plt.figure(1)
    plt.ylabel('Anzahl der Sequenzpaare')
    plt.xlabel('Größe der Kodierung in Bit')

    #maxkey = max(max(counter1.keys()),max(counter2.keys()),max(counter3.keys()))
    #maxend = max(counter1.most_common(1)[0][1],counter2.most_common(1)[0][1],
    #           counter3.most_common(1)[0][1])
    maxkey = max(max(counter1.keys()),max(counter2.keys()))
    maxend = max(counter1.most_common(1)[0][1],counter2.most_common(1)[0][1])
    plt.axis([0, maxkey * 1.2, 0, maxend * 1.2])
  
    plt.plot(counter1.keys(), counter1.values(),'bo', 
           label="Binär")
    plt.plot(counter2.keys(), counter2.values(), 'co', 
           label="Unär")
    plt.legend(loc='upper right')

    plt.figure(2)
    plt.ylabel('Anzahl der Sequenzpaare')
    plt.xlabel('Größe der Kodierung in Bit')
    plt.axis([0, max(counter3.keys()) * 1.2, 
              0, counter3.most_common(1)[0][1] * 1.2])

    plt.plot(counter3.keys(), counter3.values(), 'ro', 
           label="Huffman")
    plt.legend(loc='upper right')
  

  bin_mean = float(sum(bs1))/len(bs1)
  una_mean = float(sum(bs2))/len(bs2)
  huf_mean = float(sum(bs3))/len(bs3)

  print counter1
  print counter2
  print counter3
  print "##########################"
  print counter1.keys()
  print counter2.keys()
  print counter3.keys()

  print "Bin/Una:",float(bin_mean)/una_mean
  print "Una/Bin:",float(una_mean)/bin_mean

  print "Bin/Huf:",float(bin_mean)/huf_mean
  print "Huf/Bin:",float(huf_mean)/bin_mean

  print "Huf/Una:",float(huf_mean)/una_mean
  print "Una/Huf:",float(una_mean)/huf_mean


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

  """
  y = list(range(counter1.most_common(1)[0][1]))
  x = [bin_mean] * counter1.most_common(1)[0][1]

  print "Binary Mean:", bin_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  print "Figure 1 complete.."

  #plt.figure(2)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter2.keys())*1.2, 0, counter2.most_common(1)[0][1] * 1.2])
  plt.plot(counter2.keys(), counter2.values(), 'bo', 
           label="Unary Coding")
  una_mean = float(sum(bs2))/len(bs2)
  y = list(range(counter2.most_common(1)[0][1]))
  x = [una_mean] * counter2.most_common(1)[0][1]

  print "Unary Mean:", una_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  print "Figure 2 complete.."

  #plt.figure(3)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter3.keys())*1.2, 0, counter3.most_common(1)[0][1] * 1.2])
  plt.plot(counter3.keys(), counter3.values(), 'bo', 
           label="Huffman Coding")
  huf_mean = float(sum(bs3))/len(bs3)
  y = list(range(counter3.most_common(1)[0][1]))
  x = [huf_mean] * counter3.most_common(1)[0][1]

  print "Huffman Mean:", huf_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  print "Figure 3 complete.."

  print counter1
  print counter2
  print counter3
  print "##########################"
  print counter1.keys()
  print counter2.keys()
  print counter3.keys()

  print "Bin/Una:",float(bin_mean)/una_mean
  print "Una/Bin:",float(una_mean)/bin_mean
  print "Bin/Huf:",float(bin_mean)/huf_mean
  print "Huf/Bin:",float(huf_mean)/bin_mean
  print "Huf/Una:",float(huf_mean)/una_mean
  print "Una/Huf:",float(una_mean)/huf_mean

  print "Calculation complete.\nClock time: %.2f seconds." % (time.clock() - t)
  plt.show()
  """

def multiple_delta(a, b, c, t):
  counter1 = collections.Counter(bucket(a,3))
  counter2 = collections.Counter(bucket(b,3))
  counter3 = collections.Counter(bucket(c,3))

  mean1 = mean2 = mean3 = 0

  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter2.keys())*1.2, 0, counter2.most_common(1)[0][1] * 1.2])
  plt.plot(counter1.keys(), counter1.values(),'bo', 
           label="Delta 50")
  mean1 = float(sum(a))/len(a)
  y = list(range(counter1.most_common(1)[0][1]))
  x = [mean1] * counter1.most_common(1)[0][1]

  print "Mean1:", mean1
  plt.plot(x,y, 'b--',label='Mean1')

  plt.plot(counter2.keys(), counter2.values(),'go', 
           label="Delta 100")
  mean2 = float(sum(b))/len(b)
  y = list(range(counter2.most_common(1)[0][1]))
  x = [mean2] * counter2.most_common(1)[0][1]

  print "Mean2:", mean2
  plt.plot(x,y, 'g--',label='Mean2')

  plt.plot(counter3.keys(), counter3.values(),'ro', 
           label="Delta 200")
  mean3 = float(sum(c))/len(c)
  y = list(range(counter3.most_common(1)[0][1]))
  x = [mean3] * counter3.most_common(1)[0][1]

  print "Mean3:", mean3
  plt.plot(x,y, 'r--',label='Mean3')

  plt.legend(loc='upper right')

  print "Unterschied 1/2:", float(mean1)/mean2
  print "Unterschied 2/1:", float(mean2)/mean1
  print "Unterschied 1/3:", float(mean1)/mean3
  print "Unterschied 2/3:", float(mean2)/mean3
  print "Unterschied 3/1:", float(mean3)/mean1
  print "Unterschied 3/2:", float(mean3)/mean2

  print "Calculation complete.\nClock time: %.2f seconds." % (time.clock() - t)
  plt.show()

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
  #print "Symbol\tWeight\tHuffman Code"
  for p in huff:
    bits += len(p[1])
    it = (p[0],p[1],len(p[1]))
    canon.append(it)
    #print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
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

  keys = c.keys()
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

def megaplot(bs1, bs2, bs3, bs4, bs5, bs6, t):

  counter1 = collections.Counter(bucket(bs1,10))
  counter2 = collections.Counter(bucket(bs2,10))
  counter3 = collections.Counter(bucket(bs3,10))

  cig_bin_mean = cig_una_mean = cig_huf_mean = 0

  plt.figure(1)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter1.keys())*1.2, 0, counter1.most_common(1)[0][1] * 1.2])
  plt.plot(counter1.keys(), counter1.values(),'bo', 
           label="Binary Coding")
  cig_bin_mean = float(sum(bs1))/len(bs1)
  y = list(range(counter1.most_common(1)[0][1]))
  x = [cig_bin_mean] * counter1.most_common(1)[0][1]

  print "Cigar Binary Mean:", cig_bin_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')

  plt.figure(2)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter2.keys())*1.2, 0, counter2.most_common(1)[0][1] * 1.2])
  plt.plot(counter2.keys(), counter2.values(), 'bo', 
           label="Unary Coding")
  cig_una_mean = float(sum(bs2))/len(bs2)
  y = list(range(counter2.most_common(1)[0][1]))
  x = [cig_una_mean] * counter2.most_common(1)[0][1]

  print "Cigar Unary Mean:", cig_una_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')

  plt.figure(3)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter3.keys())*1.2, 0, counter3.most_common(1)[0][1] * 1.2])
  plt.plot(counter3.keys(), counter3.values(), 'bo', 
           label="Huffman Coding")
  cig_huf_mean = float(sum(bs3))/len(bs3)
  y = list(range(counter3.most_common(1)[0][1]))
  x = [cig_huf_mean] * counter3.most_common(1)[0][1]

  print "Cigar Huffman Mean:", cig_huf_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')



  counter4 = collections.Counter(bucket(bs4,3))
  counter5 = collections.Counter(bucket(bs5,3))
  counter6 = collections.Counter(bucket(bs6,3))

  diff_bin_mean = diff_una_mean = diff_huf_mean = 0

  plt.figure(4)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter4.keys())*1.2, 0, counter4.most_common(1)[0][1] * 1.2])
  plt.plot(counter4.keys(), counter4.values(),'bo', 
           label="Binary Coding")
  diff_bin_mean = float(sum(bs4))/len(bs4)
  y = list(range(counter4.most_common(1)[0][1]))
  x = [diff_bin_mean] * counter4.most_common(1)[0][1]

  print "Diff Binary Mean:", diff_bin_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')

  plt.figure(5)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter5.keys())*1.2, 0, counter5.most_common(1)[0][1] * 1.2])
  plt.plot(counter5.keys(), counter5.values(), 'bo', 
           label="Unary Coding")
  diff_una_mean = float(sum(bs5))/len(bs5)
  y = list(range(counter5.most_common(1)[0][1]))
  x = [diff_una_mean] * counter5.most_common(1)[0][1]

  print "Diff Unary Mean:", diff_una_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')

  plt.figure(6)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter6.keys())*1.2, 0, counter6.most_common(1)[0][1] * 1.2])
  plt.plot(counter6.keys(), counter6.values(), 'bo', 
           label="Huffman Coding")
  diff_huf_mean = float(sum(bs6))/len(bs6)
  y = list(range(counter6.most_common(1)[0][1]))
  x = [diff_huf_mean] * counter6.most_common(1)[0][1]

  print "Diff Huffman Mean:", diff_huf_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')


  print "##########################"
  print counter1
  print counter2
  print counter3
  print "##########################"
  print counter1.keys()
  print counter2.keys()
  print counter3.keys()
  print "##########################"

  print "Cig Bin/Una:",float(cig_bin_mean)/cig_una_mean
  print "Cig Una/Bin:",float(cig_una_mean)/cig_bin_mean
  print "Cig Bin/Huf:",float(cig_bin_mean)/cig_huf_mean
  print "Cig Huf/Bin:",float(cig_huf_mean)/cig_bin_mean
  print "Cig Huf/Una:",float(cig_huf_mean)/cig_una_mean
  print "Cig Una/Huf:",float(cig_una_mean)/cig_huf_mean


  print "##########################"
  print counter4
  print counter5
  print counter6
  print "##########################"
  print counter4.keys()
  print counter5.keys()
  print counter6.keys()
  print "##########################"

  print "Diff Bin/Una:",float(diff_bin_mean)/diff_una_mean
  print "Diff Una/Bin:",float(diff_una_mean)/diff_bin_mean
  print "Diff Bin/Huf:",float(diff_bin_mean)/diff_huf_mean
  print "Diff Huf/Bin:",float(diff_huf_mean)/diff_bin_mean
  print "Diff Huf/Una:",float(diff_huf_mean)/diff_una_mean
  print "Diff Una/Huf:",float(diff_una_mean)/diff_huf_mean

  print "Unterschied Bin:", float(diff_bin_mean)/cig_bin_mean
  print "Unterschied Una:", float(diff_una_mean)/cig_una_mean
  print "Unterschied Huf:", float(diff_huf_mean)/cig_huf_mean

  print "Calculation complete.\nClock time: %.2f seconds." % (time.clock() - t)
  plt.show()

def main():

  t = time.clock()
  print "Läuft..."
  amount = 100
  length = 1000
  delta = 100
  err_rate = 0.15
  
  # CIGAR

  ciglist, delta = create_cigars_or_tp(amount,length,err_rate,"acgt",delta, "cigar")

  print "LÄNGE:", len(ciglist)

  bs1 = kodierung("cigar","binary", delta, ciglist, [])
  bs2 = kodierung("cigar","unary", delta, ciglist, [])
  bs3 = kodierung("cigar","huffman", delta, ciglist, [])

  multiplot(bs1, bs2, bs3, t, 10, "cigar")
  """
  # Differenzen

  TP, delta = create_cigars_or_tp(amount,length,err_rate,"acgt",delta, "tracepoint")

  print "LÄNGE:", len(TP)

  bs1 = kodierung("tracepoint","binary", delta, [], TP)
  bs2 = kodierung("tracepoint","unary", delta, [], TP)
  bs3 = kodierung("tracepoint","huffman", delta, [], TP)

  multiplot(bs1, bs2, bs3, t, 3, "tracepoint")
  """
  """
  # Entropy
  entropy(100,1000,0.15,"acgt",100)
  """












  """ 
  bs1 = calc_bits("cigar","binary",10000,1000,0.15,"acgt",100,"cig-bin-10-1000-d100")
  bs2 = calc_bits("cigar","unary",10000,1000,0.15,"acgt",100,"cig-una-10-1000-d100")
  bs3 = calc_bits("cigar","huffman",10000,1000,0.15,"acgt",100,"cig-huf-10-1000-d100")
  
  #multiplot(bs1, bs2, bs3, t, "cig")

  thislist1 = []
  for key,value in counter11.items():
    print key, value
    for i in range(0,value):
      thislist1.append(key)
  thislist2 = []
  for key,value in counter22.items():
    print key, value
    for i in range(0,value):
      thislist2.append(key)
  thislist3 = []
  for key,value in counter33.items():
    print key, value
    for i in range(0,value):
      thislist3.append(key)

  multiplot(thislist1, thislist2, thislist3, t, "cig")
   
  bs4 = calc_bits("cigar","binary",10000,1000,0.15,"acgt",100,
            "diff-bin-10-1000-d100")

  bs5 = calc_bits("cigar","unary",10000,10000,0.15,"acgt",100,
            "diff-una-10-1000-d100")

  bs6 = calc_bits("cigar","huffman",10000,10000,0.15,"acgt",100,
            "diff-huf-10-1000-d100")
  multiplot(bs4, bs5, bs6, t, "cig")

  # multiple delta values
  del50 = calc_bits("tracepoint","unary",1000,1000,0.15,"acgt",50,
                    "diff-una-d50-1000-1000")
  del100 = calc_bits("tracepoint","unary",1000,1000,0.15,"acgt",100,
                    "diff-una-d100-1000-1000")
  del200 = calc_bits("tracepoint","unary",1000,1000,0.15,"acgt",200,
                    "diff-una-d200-1000-1000")
  multiple_delta(del50, del100, del200, t)
  """ 
  #megaplot(bs1, bs2, bs3, bs4, bs5, bs6, t)
  
  #entropy(10000,1000,0.15,"acgt",100)


if __name__ == "__main__":
  main()

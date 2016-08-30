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
#################

counter11 = ({122: 670, 366: 670, 126: 667, 378: 667, 134: 666, 469: 666, 130: 652, 455: 652, 138: 623, 483: 623, 118: 594, 354: 594, 114: 575, 342: 575, 142: 572, 497: 572, 110: 481, 330: 481, 525: 442, 150: 442, 146: 440, 511: 440, 106: 431, 318: 431, 539: 383, 154: 383, 553: 351, 158: 351, 102: 297, 306: 297, 567: 274, 162: 274, 98: 251, 294: 251, 581: 217, 166: 217, 94: 177, 282: 177, 595: 139, 170: 139, 609: 128, 174: 128, 90: 114, 270: 114, 86: 85, 258: 85, 623: 74, 178: 74, 82: 58, 246: 58, 637: 55, 182: 55, 78: 42, 234: 42, 128: 34, 651: 34, 186: 34, 384: 34, 124: 31, 372: 31, 136: 30, 476: 30, 665: 26, 190: 26, 140: 25, 144: 25, 210: 25, 490: 25, 504: 25, 120: 23, 360: 23, 70: 20, 222: 20, 532: 19, 74: 19, 132: 19, 152: 19, 462: 19, 112: 18, 336: 18, 108: 17, 160: 17, 324: 17, 560: 16, 116: 16, 348: 16, 546: 14, 156: 14, 198: 13, 588: 12, 679: 12, 168: 12, 194: 12, 574: 11, 616: 11, 164: 11, 176: 11, 100: 10, 300: 10, 518: 9, 96: 9, 148: 9, 288: 9, 92: 8, 693: 8, 707: 8, 202: 8, 276: 8, 206: 7, 721: 7, 602: 6, 172: 6, 66: 5, 88: 5, 735: 5, 264: 5, 104: 4, 312: 4, 180: 3, 630: 3, 214: 3, 218: 3, 749: 3, 763: 3, 686: 2, 80: 2, 84: 2, 155: 2, 192: 2, 196: 2, 240: 2, 252: 2, 62: 2, 672: 2, 58: 1, 64: 1, 72: 1, 145: 1, 184: 1, 700: 1, 200: 1, 216: 1, 254: 1, 777: 1, 644: 1, 889: 1})
counter22 = ({110: 2041, 120: 1833, 100: 1601, 130: 1303, 90: 1184, 140: 752, 80: 546, 410: 465, 400: 457, 420: 450, 450: 433, 440: 431, 390: 429, 430: 425, 460: 405, 380: 403, 370: 396, 360: 388, 470: 341, 480: 334, 340: 327, 150: 321, 350: 311, 330: 295, 490: 288, 310: 266, 500: 266, 320: 246, 510: 226, 520: 219, 300: 218, 530: 203, 70: 195, 290: 182, 540: 163, 550: 161, 280: 154, 270: 123, 160: 121, 560: 117, 260: 106, 250: 86, 570: 85, 580: 79, 240: 75, 590: 73, 230: 52, 170: 50, 600: 50, 610: 42, 60: 41, 620: 36, 220: 34, 180: 27, 630: 23, 210: 18, 200: 17, 640: 16, 650: 12, 660: 12, 190: 12, 680: 9, 670: 7, 50: 7, 690: 4, 700: 3, 720: 2, 770: 1, 750: 1, 760: 1})
counter33 = ({160: 2556, 170: 1564, 180: 1125, 150: 1052, 190: 865, 140: 830, 130: 615, 200: 603, 120: 381, 210: 207, 110: 64, 100: 51, 220: 50, 90: 24, 230: 10, 80: 1, 240: 1, 250: 1})


#################
def calc_bits(method, mode, amount, random_length, error_rate, alphabet, delta,
              filename):
  start_seq1 = start_seq2 = 0
  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  bit_sum = []

  cig_count_sum = []
  for i in range(0, len(random_seq_list), 2):
    end_seq1 = len(random_seq_list[i])
    end_seq2 = len(random_seq_list[i+1]) 
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                    start_seq1, end_seq1, start_seq2, end_seq2)

    cigar = aln.calc_cigar(aln.seq1, aln.seq2)
    # huffman coding for cigar
    if method == "cigar":
      if mode == "huffman":
        for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
          if cig_symbol == 'M':
            cig_count_sum.append(1)
          else:
            cig_count_sum.append(2)
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
      tp_aln = TracePoint.TracePointAlignment(aln.seq1, aln.seq2, start_seq1, 
                                            end_seq1, start_seq2, end_seq2, 
                                            delta, cigar)

      # store differences between Trace Points and Delta Value in List
      TP = [delta, tp_aln.tp[0]]
      for j in range(1,len(tp_aln.tp)):
        TP.append(tp_aln.tp[j]-tp_aln.tp[j-1])
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
    bucket.append(int(base * round(float(element)/base)))

  return bucket

def multiplot(bs1, bs2, bs3, t, method):
  counter1 = collections.Counter(bucket(bs1))
  counter2 = collections.Counter(bucket(bs2))
  counter3 = collections.Counter(bucket(bs3))

  bin_mean = una_mean = huf_mean = 0

  plt.figure(1)
  plt.ylabel('Anzahl der Sequenzpaare')
  plt.xlabel('Größe der Kodierung in Bit')
  plt.axis([0, max(counter1.keys())*1.2, 0, counter1.most_common(1)[0][1] * 1.2])
  plt.plot(counter1.keys(), counter1.values(),'bo', 
           label="Binary Coding")
  bin_mean = float(sum(bs1))/len(bs1)
  y = list(range(counter1.most_common(1)[0][1]))
  x = [bin_mean] * counter1.most_common(1)[0][1]

  print "Binary Mean:", bin_mean
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/%s-bin-1000-1000-d100" % method)

  plt.figure(2)
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
  plt.savefig("runs/%s-una-1000-1000-d100" % method)

  plt.figure(3)
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
  plt.savefig("runs/%s-huf-1000-1000-d100" % method)

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
  bits = avg = 0
  # print "Symbol\tWeight\tHuffman Code"
  for p in huff:
    bits += len(p[1])
    # print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
  return bits 

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

  cig_counter = collections.Counter(cig_ent_sum)
  diff_counter = collections.Counter(diff_ent_sum)

  cig_mean = diff_mean = 0

  plt.figure(1)
  plt.ylabel('Anzahl Sequenzpaare')
  plt.xlabel('Entropie in Bit')
  plt.axis([0, max(cig_counter.keys())*1.2, 0, 
            cig_counter.most_common(1)[0][1] * 1.2])
  cigkeylist = [int(i) for i in cig_counter.keys()]
  diffkeylist = [int(i) for i in diff_counter.keys()]
  plt.plot(cigkeylist, cig_counter.values(), 
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
  plt.plot(diffkeylist, diff_counter.values(), 
           'ro', label="Differences Entropy")
  diff_mean = float(sum(diff_ent_sum))/len(diff_ent_sum)
  y = list(range(diff_counter.most_common(1)[0][1]))
  x = [diff_mean] * diff_counter.most_common(1)[0][1]
  print "Diff Entropy Mean:", diff_mean
  plt.plot(x,y, 'r--',label='Mean')
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
   
  bs1 = calc_bits("cigar","binary",10000,1000,0.15,"acgt",100,"cig-bin-10-1000-d100")
  bs2 = calc_bits("cigar","unary",10000,1000,0.15,"acgt",100,"cig-una-10-1000-d100")
  bs3 = calc_bits("cigar","huffman",10000,1000,0.15,"acgt",100,"cig-huf-10-1000-d100")
  
  #entropy(100,1000,0.15,"acgt",100)
  #multiplot(bs1, bs2, bs3, t, "cig")

  """
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
  """
  bs4 = calc_bits("tracepoint","binary",10000,1000,0.15,"acgt",100,
            "diff-bin-10-1000-d100")

  bs5 = calc_bits("tracepoint","unary",10000,1000,0.15,"acgt",100,
            "diff-una-10-1000-d100")

  bs6 = calc_bits("tracepoint","huffman",10000,1000,0.15,"acgt",100,
            "diff-huf-10-1000-d100")
  #multiplot(bs4, bs5, bs6, t, "diff")
  
  megaplot(bs1, bs2, bs3, bs4, bs5, bs6, t)

if __name__ == "__main__":
  main()

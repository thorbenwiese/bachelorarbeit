import tp_calc
import Alignment
import TracePoint
import Cigar_Pattern
import math
import matplotlib.pyplot as plt
import collections
import argparse
import time

from heapq import heappush, heappop, heapify
from collections import defaultdict

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

def multiplot(bs1, bs2, bs3, t, method):
  counter1 = collections.Counter(bs1)
  counter2 = collections.Counter(bs2)
  counter3 = collections.Counter(bs3)

  plt.figure(1)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Anzahl Bits')
  plt.axis([0, max(counter1.keys())*1.2, 0, counter1.most_common(1)[0][1] * 1.2])
  plt.plot(counter1.keys(), counter1.values(),'bo', label="Binary Coding")
  y = list(range(counter1.most_common(1)[0][1]))
  x = [float(sum(bs1))/len(bs1)] * counter1.most_common(1)[0][1]
  print "Binary Mean:",float(sum(bs1))/len(bs1)
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/%s-bin-1000-1000-d100" % method)

  plt.figure(2)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Anzahl Bits')
  plt.axis([0, max(counter2.keys())*1.2, 0, counter2.most_common(1)[0][1] * 1.2])
  plt.plot(counter2.keys(), counter2.values(), 'bo', label="Unary Coding")
  y = list(range(counter2.most_common(1)[0][1]))
  x = [float(sum(bs2))/len(bs2)] * counter2.most_common(1)[0][1]
  print "Unary Mean:",float(sum(bs2))/len(bs2)
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/%s-una-1000-1000-d100" % method)

  plt.figure(3)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Anzahl Bits')
  plt.axis([0, max(counter3.keys())*1.2, 0, counter3.most_common(1)[0][1] * 1.2])
  plt.plot(counter3.keys(), counter3.values(), 'bo', label="Huffman Coding")
  y = list(range(counter3.most_common(1)[0][1]))
  x = [float(sum(bs3))/len(bs3)] * counter3.most_common(1)[0][1]
  print "Huffman Mean:",float(sum(bs3))/len(bs3)
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/%s-huf-1000-1000-d100" % method)


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
      #if cig_symbol == 'M':
      #  mcount += 1
      #elif cig_symbol == 'D':
      #  dcount += 1
      #elif cig_symbol == 'I':
      #  icount += 1
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
    #cig_ent_sum.append(-ent1 + (-ent2))

    # store differences between Trace Points and Delta Value in List
    TP = [delta, tp_aln.tp[0]]
    for j in range(1,len(tp_aln.tp)):
      TP.append(tp_aln.tp[j]-tp_aln.tp[j-1])

    counter = collections.Counter(TP)
    ent = 0
    for value in counter.values():
      px = float(value)/sum(counter.values())
      ent += px * math.log(px,2)
    diff_ent_sum.append(-ent*len(TP))

  cig_counter = collections.Counter(cig_ent_sum)
  diff_counter = collections.Counter(diff_ent_sum)

  plt.figure(1)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Entropie')
  plt.axis([0, max(cig_counter.keys())*1.2, 0, 
            cig_counter.most_common(1)[0][1] * 1.2])
  plt.plot(cig_counter.keys(), cig_counter.values(), 
           'bo', label="CIGAR Entropy")
  y = list(range(cig_counter.most_common(1)[0][1]))
  x = [float(sum(cig_ent_sum))/len(cig_ent_sum)] * cig_counter.most_common(1)[0][1]
  print "Cigar Entropy Mean:",float(sum(cig_ent_sum))/len(cig_ent_sum)
  plt.plot(x,y, 'b--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/cig_entropy")

  plt.figure(2)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Entropie')
  plt.axis([0, max(diff_counter.keys())*1.2, 0, 
            diff_counter.most_common(1)[0][1] * 1.2])
  plt.plot(diff_counter.keys(), diff_counter.values(), 
           'ro', label="Differences Entropy")
  y = list(range(diff_counter.most_common(1)[0][1]))
  x = [float(sum(diff_ent_sum))/len(diff_ent_sum)] * diff_counter.most_common(1)[0][1]
  print "Diff Entropy Mean:",float(sum(diff_ent_sum))/len(diff_ent_sum)
  plt.plot(x,y, 'r--',label='Mean')
  plt.legend(loc='upper right')
  plt.savefig("runs/diff_entropy")

  plt.show()
def main():

  t = time.clock()
  """ 
  bs1 = calc_bits("cigar","binary",10000,1000,0.15,"acgt",100,"cig-bin-10000-1000-d100")
  bs2 = calc_bits("cigar","unary",10000,1000,0.15,"acgt",100,"cig-una-10000-1000-d100")
  bs3 = calc_bits("cigar","huffman",10000,1000,0.15,"acgt",100,"cig-huf-10000-1000-d100")
  #entropy(1000,1000,0.15,"acgt",100)
  multiplot(bs1, bs2, bs3, t, "cig")
  """
  bs4 = calc_bits("tracepoint","binary",10000,1000,0.15,"acgt",100,
            "diff-bin-10000-1000-d100")

  bs5 = calc_bits("tracepoint","unary",10000,1000,0.15,"acgt",100,
            "diff-una-10000-1000-d100")

  bs6 = calc_bits("tracepoint","huffman",10000,1000,0.15,"acgt",100,
            "diff-huf-10000-1000-d100")
  multiplot(bs4, bs5, bs6, t, "diff")
  
if __name__ == "__main__":
  main()

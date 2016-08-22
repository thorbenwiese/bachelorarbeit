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
        # bit_sum.append(huffman(cigar))
      elif mode == "binary":
        # naiive binary coding for cigar
        # bit_sum.append(int(math.ceil(math.log(len(cigar),2)))*len(cigar))
        for cig_count, cig_symbol in Cigar_Pattern.parse_cigar(cigar):
          cig_count_sum.append(cig_count)
        bit_sum.append(int(math.ceil(math.log(len(cig_count_sum),2)))*len(
                           cig_count_sum))
        bit_sum.append(2 * len(cig_count_sum))
        cig_count_sum = []
      else:
        # unary coding
        # bit_sum.append(unary(cigar))
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
  #plot(bit_sum, mode, filename)
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

def plot(bit_sum, mode, filename):
  counter = collections.Counter(bit_sum)
  if mode == "binary":
    mlabel = "Binary Coding"
  elif mode == "unary":
    mlabel = "Unary Coding"
  elif mode == "huffman":
    mlabel = "Huffman Coding"
  plt.plot(counter.keys(), counter.values(), 'bo', label=mlabel)
  y = list(range(counter.most_common(1)[0][1]))
  x = [float(sum(bit_sum))/len(bit_sum)] * counter.most_common(1)[0][1]
  print "Mean:", float(sum(bit_sum))/len(bit_sum)
  plt.plot(x,y, 'i--',label='Mean')
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Anzahl Bits')
  plt.axis([0, max(counter.keys())*1.2, 0, counter.most_common(1)[0][1] * 1.2])
  plt.legend(loc='upper right')
  plt.savefig("runs/%s" % filename)
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

def main():

  """
  parser = argparse.ArgumentParser()
  parser.add_argument("-t", "--tracepoint", help="Test with TPs", 
                      action="store_true")
  parser.add_argument("-c", "--cigar", help="Test with CIGAR-String",
                      action="store_true")
  parser.add_argument("-hu", "--huffman", help="Test with Huffman coding",
                      action="store_true")
  parser.add_argument("-b", "--binary", help="Test with binary coding",
                      action="store_true")
  parser.add_argument("-u", "--unary", help="Test with unary coding",
                      action="store_true")
  parser.add_argument("-r", "--random",
    help="Random sequences with <Amount> <Length> <Error Rate> <Alphabet>",
    nargs=4)
  parser.add_argument("-n", "--name", help="File Name")

  parser.add_argument("-d", "--delta", help="Delta Value",
                      type=int)

  args = parser.parse_args()

  method = mode = filename = ""
  delta = args.delta
  if args.tracepoint:
    method = "tracepoint"
  else:
    method = "cigar"
  if args.huffman:
    mode = "huffman"
  elif args.binary:
    mode = "binary"
  else:
    mode = "unary"

  filename = args.name
  
  calc_bits(method, mode, int(args.random[0]), int(args.random[1]), 
            float(args.random[2]), args.random[3], delta, filename)
  """
  t = time.clock()
  """
  bs1 = calc_bits("cigar","binary",100,1000,0.15,"acgt",100,"cig-bin-1000-1000-d100")
  bs2 = calc_bits("cigar","unary",100,1000,0.15,"acgt",100,"cig-una-1000-1000-d100")
  bs3 = calc_bits("cigar","huffman",100,1000,0.15,"acgt",100,"cig-huf-1000-1000-d100")

  multiplot(bs1, bs2, bs3, t)
  """
  bs4 = calc_bits("tracepoint","binary",1000,1000,0.15,"acgt",100,
            "diff-bin-1000-1000-d100")

  bs5 = calc_bits("tracepoint","unary",1000,1000,0.15,"acgt",100,
            "diff-una-1000-1000-d100")

  bs6 = calc_bits("tracepoint","huffman",1000,1000,0.15,"acgt",100,
            "diff-huf-1000-1000-d100")
  multiplot(bs4, bs5, bs6, t, "diff")

if __name__ == "__main__":
  main()

import tp_calc
import Alignment
import TracePoint
import Cigar_Pattern
import math
import matplotlib.pyplot as plt
import collections
import argparse

from heapq import heappush, heappop, heapify
from collections import defaultdict

def calc_bits(method, mode, amount, random_length, error_rate, alphabet, delta,
              filename):
  start_seq1 = start_seq2 = 0
  random_seq_list = tp_calc.random_sequences(amount,random_length,error_rate,
                    alphabet)

  bit_sum = []
  for i in range(0, len(random_seq_list), 2):
    end_seq1 = len(random_seq_list[i])
    end_seq2 = len(random_seq_list[i+1]) 
    aln = Alignment.Alignment(random_seq_list[i], random_seq_list[i + 1], 
                    start_seq1, end_seq1, start_seq2, end_seq2)

    cigar = aln.calc_cigar(aln.seq1, aln.seq2)
    # huffman coding for cigar
    if method == "cigar":
      if mode == "huffman":
        bit_sum.append(huffman(cigar))
      elif mode == "binary":
        # naiive binary coding for cigar
        bit_sum.append(int(math.ceil(math.log(len(cigar),2)))*len(cigar))
      else:
        # unary coding
        bit_sum.append(unary(cigar))
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
  plot(bit_sum, mode, filename)

def plot(bit_sum, mode, filename):
  counter = collections.Counter(bit_sum)
  if mode == "binary":
    mlabel = "Binary Coding"
  elif mode == "unary":
    mlabel = "Unary Coding"
  elif mode == "huffman":
    mlabel = "Huffman Coding"
  plt.plot(counter.keys(), counter.values(), 'bo', label=mlabel)
  plt.plot(sum(bit_sum)/len(bit_sum), 'r--')
  y = list(range(counter.most_common(1)[0][1]))
  x = [float(sum(bit_sum))/len(bit_sum)] * counter.most_common(1)[0][1]
  plt.plot(x,y, 'g--',label='Mean')
  plt.savefig("runs/%s" % filename)
  plt.ylabel('Anzahl Ergebnisse')
  plt.xlabel('Anzahl Bits')
  plt.axis([0, max(counter.keys())*1.2, 0, counter.most_common(1)[0][1] * 1.2])
  plt.legend(loc='upper right')
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

if __name__ == "__main__":
  main()

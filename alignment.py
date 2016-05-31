#!/bin/env python
# -*- coding: utf-8 -*-

import math
import re
import sys
import os
import subprocess
import argparse
import pysam
import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

# Alignment Struktur
class Alignment(object):
    def __init__(self, s1, mid, s2, score):
        self.s1 = s1
        self.mid = mid
        self.s2 = s2
        self.score = score


# Trace Point Struktur
class TracePoint(object):
    def __init__(self, uvalue, vvalue):
        self.uvalue = uvalue
        self.vvalue = vvalue


# Character / InDel Count Struktur
class InDel(object):
    def __init__(self, letter_count, indel_count):
        self.letter_count = letter_count
        self.indel_count = indel_count


# Aus Cigar und Sequenzen Alignment bauen
def cig_to_align(seq1, seq2, cig):
    # 3 Zeilen für Alignment
    seq1_align = ""
    middle = ""
    seq2_align = ""
    # tmp 1 und 2 für die Erweiterung der Sequenzen mit '-' für InDels
    tmp1 = 0
    tmp2 = 0
    # Anzahl der Edit-Operationen im Cigar-String
    cig_count = 0
    # neues Pattern für Cigar-Strings im Format: Zahl + 1 Buchstabe aus {M,I,D,N,S,H,P}
    cigar_pattern = re.compile(r"\d+[MIDNSHP=j]{1}")
    # Suche in cig nach Pattern
    for j in cigar_pattern.findall(cig):
        tmp1 += cig_count
        tmp2 += cig_count
        cig_count = int(j[:-1])
        cig_symbol = j[-1]
        if cig_symbol == 'M':
            seq1_align += str(seq1[tmp1:tmp1 + cig_count])
            middle += str("|" * cig_count)
            seq2_align += str(seq2[tmp2:tmp2 + cig_count])
        if cig_symbol == 'I':
            seq1_align += str(seq1[tmp1:tmp1 + cig_count])
            middle += str(" " * cig_count)
            seq2_align += str("-" * cig_count)
            tmp2 -= cig_count
        if cig_symbol == 'D':
            seq1_align += str("-" * cig_count)
            middle += str(" " * cig_count)
            seq2_align += str(seq2[tmp2:tmp2 + cig_count])
            tmp1 -= cig_count

    # Speicher Alignment in Alignment Struktur
    alignment = Alignment(seq1_align, middle, seq2_align,0)

    print "Gesamtalignment:\n"
    print numbers(alignment.s1)
    print alignment.s1
    print alignment.mid
    print alignment.s2
    print numbers(alignment.s2)
    print ""
    return alignment

# Berechnung des Alignments 
def calculate_alignment(sequence1, sequence2):

    # BLOSUM62 Matrix
    matrix = matlist.blosum62
    # Gap Open Penalty
    gap_open = -10
    # Gap Extend Penalty
    gap_extend = -0.5
     
    alns = pairwise2.align.globalds(sequence1.upper(), sequence2.upper(), matrix, gap_open, gap_extend)
     
    top_aln = alns[0]
    aln1, aln2, score, begin, end = top_aln

    alignment = Alignment(aln1.lower(),"",aln2.lower(),score)

    return alignment


# Entfernung von '-'
def align_to_seq(align):
    seq = align.replace("-", "")
    return seq


# TPs
def calculate_tps(seq):
    v_tp.append(seq[-1])
    return


# Anzahl an Buchstaben in Sequenz
def count_letters(seq):
    letter_count = 0
    indel_count = 0
    for k in seq:
        if k.isalpha():
            letter_count += 1
        elif k == '-':
            indel_count += 1
    counters = InDel(letter_count, indel_count)
    return counters

# 0en und 5en für das Alignment
def numbers(seq):

    numbers = ""
    for i in range(0, len(seq), 10):
            numbers += "0" + " " * 4 + "5" + " " * 4
    return numbers

# zufälliger Sequenz-Generator
def random_sequences(amount, random_length, error_rate):
    seq = ""
    random_seqs = []
    for i in range(0, amount):
        seq = ''.join(random.choice("acgt") for i in range(random_length))
        s = list(seq)

        position = 0
        for element in s:
            # random Zahl zwischen 0 und 1
            r_num = random.random()
            if r_num <= error_rate:
                s[position] = random.choice("acgt")
            position += 1
        seq_new = "".join(s)
        random_seqs.append(seq)
        random_seqs.append(seq_new)
    return random_seqs

# Berechne Striche für Alignment
def calculate_middle(seq1, seq2):
    middle = ""
    c = 0
    add = 0

    if len(seq1) < len(seq2):
        add = len(seq2) - len(seq1)
    else:
        add = len(seq1) - len(seq2)

    seq2 += " " * add
    for i in seq1:
        if i == seq2[c]:
            middle += '|'
        else:
            middle += ' '
        c += 1
    return middle


# mit den TracePoints, Sequenzen und Delta die Intervalle rekonstruieren
def rebuild_intervals(tp, seq1, seq2, delta, score, verbose):
    tp_count = int(math.ceil(float(len(seq1)) / delta)) - 1
    start_seq1 = start_seq2 = 0
    end_seq1 = start_seq1 + len(seq1) - 1
    end_seq2 = start_seq2 + len(seq2) - 1
    c = 1
    test1 = ""
    test2 = ""
    score_new = 0
    aln1 = ""
    aln2 = ""

    print "Berechnung der Intervalle anhand der Trace Points:\n"

    for i in tp.vvalue:
        if i == tp.vvalue[0]:
            test1 += str(seq1[start_seq1:delta])
            test2 += str(seq2[start_seq2:i + 1])
            aln = calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i + 1])
            aln1 += str(aln.s1)
            aln2 += str(aln.s2)
            score_new += aln.score
            if verbose:
                print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start_seq1, delta - 1, start_seq2, i)
                print aln.s1
                print calculate_middle(aln.s1,aln.s2)
                print aln.s2
                print "Score:",aln.score
                print ""
        elif i == tp.vvalue[-1]:
            test1 += str(seq1[c * delta:c * delta + delta])
            test2 += str(seq2[tp.vvalue[c - 1] + 1:i + 1])
            aln = calculate_alignment(seq1[c * delta:c * delta + delta],seq2[tp.vvalue[c - 1] + 1:i + 1])
            aln1 += str(aln.s1)
            aln2 += str(aln.s2)
            score_new += aln.score
            if verbose:
                print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (c * delta, c * delta + delta - 1, tp.vvalue[c - 1] + 1, i)
                print aln.s1
                print calculate_middle(aln.s1,aln.s2)
                print aln.s2
                print "Score:",aln.score
                print ""
            c += 1
            test1 += str(seq1[tp_count * delta:end_seq1 + 1])
            test2 += str(seq2[tp.vvalue[c - 1] + 1:end_seq2 + 1])
            aln = calculate_alignment(seq1[tp_count * delta:end_seq1 + 1],seq2[tp.vvalue[c - 1] + 1:end_seq2 + 1])
            aln1 += str(aln.s1)
            aln2 += str(aln.s2)
            score_new += aln.score
            if verbose:
                print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (tp_count * delta, end_seq1, tp.vvalue[c - 1] + 1, end_seq2)
                print aln.s1
                print calculate_middle(aln.s1,aln.s2)
                print aln.s2
                print "Score:",aln.score
                print ""

        else:
            test1 += str(seq1[c * delta:c * delta + delta])
            test2 += str(seq2[tp.vvalue[c - 1] + 1:i + 1])
            aln = calculate_alignment(seq1[c * delta:c * delta + delta],seq2[tp.vvalue[c - 1] + 1:i + 1])
            aln1 += str(aln.s1)
            aln2 += str(aln.s2)
            score_new += aln.score
            if verbose:
                print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (c * delta, c * delta + delta - 1, tp.vvalue[c - 1] + 1, i)
                print aln.s1
                print calculate_middle(aln.s1,aln.s2)
                print aln.s2
                print "Score:",aln.score
                print ""
            c += 1

    if test1 == seq1 and test2 == seq2:
        print "Die Konkatenation der Teilalignments ergibt das oben genannte Gesamtalignment!\n"
        print numbers(aln1)
        print aln1
        print calculate_middle(aln1,aln2)
        print aln2
        print numbers(aln2)
        print ""
        if score == score_new:
            print "Der Score des neuen Alignments (%.1f) ist so groß wie der des ursprünglichen Alignments (%.1f)" % (score_new, score)
        elif score < score_new:
            print "Der Score des neuen Alignments (%.1f) ist größer als der des ursprünglichen Alignments (%.1f)" % (score_new, score)
    else:
        sys.stderr.write("Falsche Rekonstruktion des Alignments!")
        if verbose:
            print seq1
            print test1
            print ""
            print seq2
            print test2
        sys.exit(1);
    return


# print Intervalle
def print_intervals(seq1, seq2, alignment, delta, verbose):
    m = 0
    v_id = 0
    start_seq1 = 0
    end_seq1 = start_seq1 + len(seq1) - 1
    start_seq2 = 0
    end_seq2 = start_seq2 + len(seq2) - 1

    # neues Alignment
    seq1_new = seq2_new = middle_new = ""

    # u_TP Liste
    u_tp = []

    # v_TP Liste
    v_tp = []

    # dynamische Berechnung von p für Intervallgröße
    p = 1
    while (p * delta) < start_seq1:
        p += 1

    # Anzahl der Intervalle
    tau = math.ceil(float(len(seq1)) / delta)
    tau = int(tau)

    # Anzahl der Trace Points
    tp_count = tau - 1

    # Intervalle
    intervals = [0] * (tp_count + 1)

    # TPs
    for q in range(0, tau):
        if q == 0:
            intervals[q] = start_seq1, p * delta - 1
        elif q == tau - 1:
            intervals[q] = (p + tau - 2) * delta, end_seq1
        else:
            intervals[q] = (p + q - 1) * delta, (p + q) * delta - 1

    # Endpunkte der Intervalle
    for i in intervals:
        u_tp.append(i[1])

    seq1_align = alignment.s1
    middle = alignment.mid
    seq2_align = alignment.s2
    score = alignment.score
    start1 = start_seq1
    start2 = start_seq2
    for n in u_tp:
        if n == u_tp[-1]:
            print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, len(seq1) - 1, start2, len(seq2) - 1)
            if verbose:
                print seq1_align[start_seq1:len(seq1_align)]
                print calculate_middle(seq1_align[start_seq1:len(seq1_align)], seq2_align[start_seq1:len(seq1_align)])
                print seq2_align[start_seq1:len(seq1_align)]
            seq1_new += str(seq1_align[start_seq1:len(seq1_align)])
            middle_new += str(middle[start_seq1:len(seq1_align)])
            seq2_new += str(seq2_align[start_seq1:len(seq1_align)])

        else:
            while count_letters(seq1_align[start_seq1:m]).letter_count != delta:
                m += 1
            indel_count = 0
            indel_count = count_letters(seq2_align[start_seq1:m]).indel_count
            v_id += indel_count
            v_tp.append(m - 1 - v_id)
            end = m - 1 - v_id
            indel_count = 0
            print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, n, start2, end)
            if verbose:
                print seq1_align[start_seq1:m]
                print calculate_middle(seq1_align[start_seq1:m], seq2_align[start_seq1:m])
                print seq2_align[start_seq1:m]
                print ""
            seq1_new += str(seq1_align[start_seq1:m])
            middle_new += str(middle[start_seq1:m])
            seq2_new += str(seq2_align[start_seq1:m])
            start_seq1 = m
            m = start_seq1
            start1 = n + 1
            start2 = start_seq1 - v_id
            end = m - 1 - v_id
    alignment_new = Alignment(seq1_new, middle_new, seq2_new,0)

    # check_alignment(alignment, alignment_new)
    tp = TracePoint(0, v_tp)
    print "\nTrace Points:", tp.vvalue, "\n"
    if check_alignment(alignment, alignment_new):
        rebuild_intervals(tp, seq1, seq2, delta, score, verbose)
    else:
        sys.stderr.write("Falsche Berechnung der Subsequenzen!")
        sys.exit(1);
    return tp


# Vergleich von Input und konkateniertem Output
def check_alignment(alignment, alignment_new):
    if alignment_new.s1 == alignment.s1 and alignment_new.s2 == alignment.s2:
        return True
    else:
        sys.stderr.write("Falsche Aufteilung der Sequenzen!")
        print ""
        print alignment_new.s1
        print alignment.s1
        print ""
        print alignment_new.s2
        print alignment.s2
        sys.exit(1);


def main(argv):
    seq1 = ""
    seq2 = ""
    cigar = ""
    delta = 0

    # Argumente
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq1", "--seq1", help="Die erste Sequenz")
    parser.add_argument("-seq2", "--seq2", help="Die zweite Sequenz")
    parser.add_argument("-d", "--delta", help="Delta", type=int)
    group1 = parser.add_mutually_exclusive_group()
    group2 = parser.add_mutually_exclusive_group()
    group1.add_argument("-c", "--cigar", help="Der CIGAR-String")
    group1.add_argument("-b", "--bam", help="Input BAM-File")
    group1.add_argument("-s", "--sam", help="Input SAM-File")
    group1.add_argument("-r", "--random", help="Zufällige Sequenzen generieren mit <Anzahl> <Länge> <Fehlerrate>",
                        nargs=3)
    group2.add_argument("-e", "--encode", help="Alignment nur mit TracePoints codieren", action="store_true")
    group2.add_argument("-x", "--extract", help="Aus TracePoints neues Alignment konstruieren", default=False, action="store_true")
    parser.add_argument("-v", "--verbose", help="Ausführlicher Output", action="store_true")
    args = parser.parse_args()

    seq1 = args.seq1
    seq2 = args.seq2
    cigar = args.cigar
    delta = args.delta
    verbose = args.verbose

    if args.verbose:
        verbose = True

    if args.seq1 and args.seq2 and args.delta and args.cigar:
        print '\nSequenz 1:', seq1
        print 'Sequenz 2:', seq2
        print 'Cigar-String:', cigar
        print 'Delta:', delta
        print ""
        alignment = cig_to_align(seq1, seq2, cigar)
        tp = print_intervals(seq1, seq2, alignment, delta, verbose)

    elif args.seq1 and args.seq2 and args.delta and not args.cigar:
        print '\nSequenz 1:', seq1
        print 'Sequenz 2:', seq2
        print 'Delta:', delta
        print ""
        alignment = calculate_alignment(seq1,seq2)

        print "Gesamtalignment:\n"
        print numbers(alignment.s1)
        print alignment.s1
        print calculate_middle(alignment.s1,alignment.s2)
        print alignment.s2
        print numbers(alignment.s2)
        print ""
        tp = print_intervals(seq1,seq2,alignment,delta,verbose)

    # TODO BAM-File noch nicht vollständig
    elif args.bam:
        bamFile = args.bam
        print "BAM-File:", bamFile
        bamFP = pysam.AlignmentFile(bamFile, "rb")
        proc = subprocess.Popen(["samtools view %s| cut -f 10 " % bamFile], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out_split = out.split('\n')
        c = 0

        for s in out_split:
            c += 1
            if c%2 == 0:
                seq2 = str(s)
            else:
                seq1 = str(s)
            alignment = cig_to_align(seq1,seq2,cigar)
            tp = print_intervals(seq1,seq2,alignment,delta)
            print "Trace Points:",tp.vvalue

        for read in bamFP:
            if( not( read.is_unmapped ) ):
                cig_line = read.cigar;
                for (cig_type,cig_length) in cig_line:
                    while z < 10:
                        try:
                            if(cig_type == 0): # match
                                cigar += "%dM" % cig_length
                            elif(cig_type == 1): # insertions
                                cigar += "%dI" % cig_length
                            elif(cig_type == 2): # deletion
                                cigar += "%dD" % cig_length
                            elif(cig_type == 3): # skip
                                cigar += "%dN" % cig_length
                            elif(cig_type == 4): # soft clipping
                                cigar += "%dS" % cig_length
                            elif(cig_type == 5): # hard clipping
                                cigar += "%dH" % cig_length
                            elif(cig_type == 6): # padding
                                cigar += "%dP" % cig_length
                            else:
                                print "Falsche Cigar-Nummer!";
                                sys.exit(1);
                            for s in out_split:
                                c += 1
                                if c%2 == 0:
                                    seq2 = str(s)
                                else:
                                    seq1 = str(s)
                                alignment = cig_to_align(seq1,seq2,cigar)
                                tp = print_intervals(seq1,seq2,alignment,delta)
                                print "Trace Points:",tp.vvalue
                        except:
                            print "Fehlerhafter Cigar-String";

    # TODO SAM-File
    elif args.sam:
        samFile = args.sam
        samFP = pysam.Samfile(bamFile, "rb")
        print "SAM-File:", samFile

    # Zufällige Sequenzen
    elif args.random:

        random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]))

        for i in range(0,int(args.random[0])*2,2):
            aln = calculate_alignment(random_seq_list[i],random_seq_list[i+1])
            if verbose:
                print aln.s1
                print calculate_middle(aln.s1, aln.s2)
                print aln.s2
                print "Score:", aln.score
                print ""
            score = aln.score

            tp = print_intervals(random_seq_list[i],random_seq_list[i+1],aln,delta,verbose)
    else:
        sys.stderr.write("Falsche Eingabe der Argumente!")
        sys.exit(1);

if __name__ == "__main__":
    main(sys.argv[1:])

#!/bin/env python
# -*- coding: utf-8 -*-

import math
import re
import sys
import os
# für Subprozess in der Konsole
import subprocess
import argparse
import pysam
import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

sys.path.append("biopython/")
sys.path.append("pysam/")

# Alignment Struktur
class TracePointAlignment(object):

    # Konstruktor
    def __init__(self, seq1, seq2, delta, score = None, tp = None, start_seq1 = 0, start_seq2 = 0):
        self.seq1 = seq1
        self.start_seq1 = start_seq1
        self.seq2 = seq2
        self.start_seq2 = start_seq2
        self.delta = delta
        self.score = score
        self.tp = tp

    # Berechnung des Alignments 
    def calculate_alignment(self, sequence1, sequence2):

        # x = No parameters.  Identical characters have score of 1, otherwise 0.
        # x = No gap penalties.
        alns = pairwise2.align.globalxx(sequence1.upper(), sequence2.upper())
        
        if len(alns) > 0:
            top_aln = alns[0]

            aln1, aln2, score, begin, end = top_aln

            alignment = TracePointAlignment(aln1.lower(),aln2.lower(),self.delta,score)
        else:
            sys.stderr.write("Es konnte kein Alignment berechnet werden.\n")
            sys.exit(1);
        
        return alignment

    # Anzahl an Buchstaben in Sequenz
    def count_indels_letters(self, seq):

        letter_count = 0
        indel_count = 0
        indel_letter_count = []
        for i in seq:
            if i.isalpha():
                letter_count += 1
            elif i == '-':
                indel_count += 1
        indel_letter_count.append(letter_count)
        indel_letter_count.append(indel_count)
        return indel_letter_count

    # Aus Cigar und Sequenzen Alignment bauen
    def cigar_to_alignment(self, seq1, seq2, cigar):
        # Zeilen für Alignment
        seq1_align = ""
        seq2_align = ""
        # tmp 1 und 2 für die Erweiterung der Sequenzen mit '-' für InDels
        tmp1 = 0
        tmp2 = 0
        # Anzahl der Edit-Operationen im Cigar-String
        cig_count = 0
        # neues Pattern für Cigar-Strings im Format: Zahl + 1 Buchstabe aus {M,I,D,N,S,H,P}
        cigar_pattern = re.compile(r"\d+[MIDNSHP=j]{1}")
        # Suche in cig nach Pattern
        for j in cigar_pattern.findall(cigar):
            tmp1 += cig_count
            tmp2 += cig_count
            cig_count = int(j[:-1])
            cig_symbol = j[-1]
            if cig_symbol == 'M':
                seq1_align += str(seq1[tmp1:tmp1 + cig_count])
                seq2_align += str(seq2[tmp2:tmp2 + cig_count])
            if cig_symbol == 'I':
                seq1_align += str(seq1[tmp1:tmp1 + cig_count])
                seq2_align += str("-" * cig_count)
                tmp2 -= cig_count
            if cig_symbol == 'D':
                seq1_align += str("-" * cig_count)
                seq2_align += str(seq2[tmp2:tmp2 + cig_count])
                tmp1 -= cig_count

        # Speicher Alignment in Alignment Struktur
        alignment = TracePointAlignment(seq1_align, seq2_align, self.delta, self.score)

        return alignment

    # 0en und 5en für das Pretty Print Alignment
    def print_sequence_positions(self, seq):

        positions = ""
        count = 0
        for i in range(0, len(seq), 5):
            if count == 0:
                positions += "0" + " " * 4
                count = 1
            elif count == 1:
                positions += "5" + " " * 4
                count = 0

        return positions


    # Ausführliche Ausgabe des Alignments mit Zahlen und Strichen
    def pretty_print_alignment(self, seq1, seq2):

        middle = ""
        count = 0
        add = 0
        pretty_print = ""

        if len(seq1) < len(seq2):
            add = len(seq2) - len(seq1)
        else:
            add = len(seq1) - len(seq2)

        seq2 += " " * add
        for i in seq1:
            if i == seq2[count]:
                middle += '|'
            else:
                middle += ' '
            count += 1

        pos1 = self.print_sequence_positions(seq1)
        pos2 = self.print_sequence_positions(seq2)
        pretty_print = "" + pos1 + "\n" + seq1 + "\n" + middle + "\n" + seq2 + "\n" + pos2
        return pretty_print

    # berechne Intervalle
    def calculate_intervals(self, alignment, verbose=False):
        m = 0
        v_id = 0
        start_seq1 = self.start_seq1
        end_seq1 = start_seq1 + len(self.seq1) - 1
        start_seq2 = self.start_seq2
        end_seq2 = start_seq2 + len(self.seq2) - 1

        # neues Alignment
        seq1_new = seq2_new = middle_new = ""

        # TPs in Sequenz 1
        tp1 = []

        # TPs in Sequenz 2
        tp2 = []

        # dynamische Berechnung der Intervallgröße
        dynamic = int(math.ceil(start_seq1/self.delta))
        if dynamic == 0:
            dynamic = 1

        # Anzahl der Intervalle
        letters_in_seq1 = self.count_indels_letters(self.seq1)[0]
        letters_in_seq2 = self.count_indels_letters(self.seq2)[0]
        # Anpassung an Länge der Sequenzen
        if letters_in_seq1 < letters_in_seq2:
            interval_count = int(math.ceil(float(letters_in_seq1) / self.delta))
        else:
            interval_count = int(math.ceil(float(letters_in_seq2) / self.delta))

        # Anzahl der Trace Points
        tp_count = interval_count - 1

        # Intervalle
        intervals = [0] * (tp_count + 1)

        # TPs
        for i in range(0, interval_count):
            # erstes Intervall
            if i == 0:
                intervals[i] = start_seq1, dynamic * self.delta - 1
            # letztes Intervall
            elif i == interval_count - 1:
                intervals[i] = (dynamic + interval_count - 2) * self.delta, end_seq1
            # restlichen Intervalle
            else:
                intervals[i] = (dynamic + i - 1) * self.delta, (dynamic + i) * self.delta - 1

        # Endpunkte der Intervalle
        for j in intervals:
            tp1.append(j[1])

        seq1_align = alignment.seq1
        seq2_align = alignment.seq2
        score = self.score
        start1 = start_seq1
        start2 = start_seq2
        for k in tp1:
            if k == tp1[-1]:
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, letters_in_seq1 - 1, start2, letters_in_seq2 - 1)#len(self.seq1/2)
                    print self.pretty_print_alignment(seq1_align[start_seq1:len(seq1_align)], seq2_align[start_seq1:len(seq1_align)])
                seq1_new += str(seq1_align[start_seq1:len(seq1_align)])
                seq2_new += str(seq2_align[start_seq1:len(seq1_align)])

            else:
                # first delta letters in seq1_align
                while self.count_indels_letters(seq1_align[start_seq1:m])[0] != self.delta:
                        # Abbruchbedingung, falls am Ende nur noch Gaps sind
                        if self.count_indels_letters(seq1_align[m:end_seq1])[0] == 0:
                            break;
                        else:
                            m += 1
                indel_count = 0
                indel_count = self.count_indels_letters(seq2_align[start_seq1:m])[1]
                v_id += indel_count
                tp2.append(m - 1 - v_id)
                end = m - 1 - v_id
                indel_count = 0
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, k, start2, end)
                    print self.pretty_print_alignment(seq1_align[start_seq1:m], seq2_align[start_seq1:m])
                    print ""
                seq1_new += str(seq1_align[start_seq1:m])
                seq2_new += str(seq2_align[start_seq1:m])
                start_seq1 = m
                m = start_seq1
                start1 = k + 1
                start2 = start_seq1 - v_id
                end = m - 1 - v_id
        TP_alignment = TracePointAlignment(seq1_new, seq2_new, self.delta, self.score, tp2)
        # Test auf Korrektheit der Sequenzen
        if self.check_alignment(alignment, TP_alignment):
            if verbose:
                print "\nKonkateniertes Alignment:"
                self.show_TracePointAlignment(TP_alignment)
            return TP_alignment


    # Ausgabe des TracePointAlignments
    def show_TracePointAlignment(self, tp_alignment):

        print self.pretty_print_alignment(tp_alignment.seq1,tp_alignment.seq2)

    # Vergleich von Input und konkateniertem Output
    def check_alignment(self, alignment1, alignment2):
        if alignment2.seq1 == alignment1.seq1 and alignment2.seq2 == alignment1.seq2:
            return True
        else:
            sys.stderr.write("Falsche Aufteilung der Sequenzen!\n")
            # für Debugging
            print ""
            print alignment1.seq1
            print alignment2.seq1
            print ""
            print alignment1.seq2
            print alignment2.seq2
            sys.exit(1);

    # mit den TracePoints, Sequenzen und Delta die Intervalle rekonstruieren
    def rebuild_intervals(self, tp_alignment, seq1, seq2, delta, score, verbose):
        tp_count = len(tp_alignment.tp)
        start_seq1 = start_seq2 = 0
        end_seq1 = start_seq1 + len(seq1) - 1
        end_seq2 = start_seq2 + len(seq2) - 1
        count = 1
        check_seq1 = check_seq2 = aln1 = aln2 = ""
        score_new = 0

        # Anzahl der Intervalle
        letters_in_seq1 = self.count_indels_letters(seq1)[0]
        letters_in_seq2 = self.count_indels_letters(seq2)[0]

        if verbose:
            print "\nBerechnung der Intervalle anhand der Trace Points:\n"

        alignment = TracePointAlignment(seq1,seq2,delta,score)

        for i in tp_alignment.tp:
            if i == tp_alignment.tp[0]:
                check_seq1 += str(seq1[start_seq1:delta])
                check_seq2 += str(seq2[start_seq2:i + 1])
                aln = alignment.calculate_alignment(seq1[start_seq1:delta],seq2[start_seq2:i + 1])
                aln1 += str(aln.seq1)
                aln2 += str(aln.seq2)
                if aln.score != None:
                    score_new += aln.score
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start_seq1, delta - 1, start_seq2, i)
                    print alignment.pretty_print_alignment(aln.seq1,aln.seq2)
                    print "Score:",aln.score
                    print ""
            elif i == tp_alignment.tp[-1]:
                check_seq1 += str(seq1[count * delta:count * delta + delta])
                check_seq2 += str(seq2[tp_alignment.tp[count - 1] + 1:i + 1])
                aln = alignment.calculate_alignment(seq1[count * delta:count * delta + delta],seq2[tp_alignment.tp[count - 1] + 1:i + 1])
                aln1 += str(aln.seq1)
                aln2 += str(aln.seq2)
                if aln.score != None:
                    score_new += aln.score
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (count * delta, count * delta + delta - 1, tp_alignment.tp[count - 1] + 1, i)
                    print alignment.pretty_print_alignment(aln.seq1,aln.seq2)
                    print "Score:",aln.score
                    print ""
                count += 1
                check_seq1 += str(seq1[tp_count * delta:end_seq1 + 1])
                check_seq2 += str(seq2[tp_alignment.tp[count - 1] + 1:end_seq2 + 1])
                aln = alignment.calculate_alignment(seq1[tp_count * delta:end_seq1 + 1],seq2[tp_alignment.tp[count - 1] + 1:end_seq2 + 1])
                aln1 += str(aln.seq1)
                aln2 += str(aln.seq2)
                if aln.score != None:
                    score_new += aln.score
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (tp_count * delta, letters_in_seq1 - 1, tp_alignment.tp[count - 1] + 1, letters_in_seq2 - 1)
                    print alignment.pretty_print_alignment(aln.seq1, aln.seq2)
                    print "Score:",aln.score
                    print ""

            else:
                check_seq1 += str(seq1[count * delta:count * delta + delta])
                check_seq2 += str(seq2[tp_alignment.tp[count - 1] + 1:i + 1])
                aln = alignment.calculate_alignment(seq1[count * delta:count * delta + delta],seq2[tp_alignment.tp[count - 1] + 1:i + 1])
                aln1 += str(aln.seq1)
                aln2 += str(aln.seq2)
                if aln.score != None:
                    score_new += aln.score
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (count * delta, count * delta + delta - 1, tp_alignment.tp[count - 1] + 1, i)
                    print alignment.pretty_print_alignment(aln.seq1, aln.seq2)
                    print "Score:",aln.score
                    print ""
                count += 1

        if check_seq1 == seq1 and check_seq2 == seq2:
            if verbose:
                print "Die Konkatenation der Teilalignments ergibt das oben genannte Gesamtalignment!\n"
                print alignment.pretty_print_alignment(aln1,aln2)
                print ""
                if score_new != None and score != None:
                    if score == score_new:
                        print "Der Score des neuen Alignments (%.1f) ist so groß wie der des ursprünglichen Alignments (%.1f)" % (score_new, score)
                    elif score < score_new:
                        print "Der Score des neuen Alignments (%.1f) ist größer als der des ursprünglichen Alignments (%.1f)" % (score_new, score)
                    else:
                        print "Der Score des neuen Alignments (%.1f) ist kleiner als der des ursprünglichen Alignments (%.1f)" % (score_new, score)
                else:
                    print "Es konnte kein neuer Score berechnet werden."
        else:
            sys.stderr.write("Falsche Rekonstruktion des Alignments!\n")
            # für Debugging
            if verbose:
                print seq1
                print check_seq1
                print ""
                print seq2
                print check_seq2
            sys.exit(1);
        return

# zufälliger Sequenz-Generator
def random_sequences(amount, random_length, error_rate, alphabet):
    seq = ""
    random_seqs = []
    for i in range(0, amount):
        seq = ''.join(random.choice(alphabet) for j in range(random_length))
        r_seq1 = list(seq)
        r_seq2 = list(seq)

        position = 0
        for element in r_seq2:
            # random Zahl zwischen 0 und 1
            r_num = random.random()
            if r_num <= error_rate:
                r_choice = random.random()
                if 0 <= r_choice < 0.1:
                    # Insertion
                    r_seq2.insert(position, '-')
                elif 0.1 < r_choice < 0.2:
                    # Deletion
                    r_seq1.insert(position, '-')
                elif 0.2 < r_choice < 0.4:
                    # a
                    r_seq2[position] = 'a'
                elif 0.4 < r_choice < 0.6:
                    # c
                    r_seq2[position] = 'c'
                elif 0.6 < r_choice < 0.8:
                    # g
                    r_seq2[position] = 'g'
                elif 0.8 < r_choice <= 1.0:
                    # t
                    r_seq2[position] = 't'
                r_seq2[position] = random.choice(alphabet)
            position += 1
        seq_new1 = "".join(r_seq1)
        seq_new2 = "".join(r_seq2)
        random_seqs.append(seq_new1)
        random_seqs.append(seq_new2)
    return random_seqs

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
    group1.add_argument("-r", "--random", help="Zufällige Sequenzen generieren mit <Anzahl> <Länge> <Fehlerrate> <Alphabet>",
                        nargs=4)
    group2.add_argument("-e", "--encode", help="Alignment nur mit TracePoints codieren", default=False, action="store_true")
    group2.add_argument("-x", "--decode", help="Aus TracePoints neues Alignment konstruieren", default=False, action="store_true")
    parser.add_argument("-v", "--verbose", help="Ausführlicher Output", default=False, action="store_true")
    args = parser.parse_args()

    seq1 = args.seq1
    seq2 = args.seq2
    cigar = args.cigar
    delta = args.delta
    verbose = args.verbose
    decode = args.decode

    # Mit CIGAR
    if args.seq1 and args.seq2 and args.delta and args.cigar:
        print '\nSequenz 1:', seq1
        print 'Sequenz 2:', seq2
        print 'Cigar-String:', cigar
        print 'Delta:', delta
        print ""
        alignment = TracePointAlignment(seq1,seq2,delta)
        alignment1 = alignment.cigar_to_alignment(alignment.seq1, alignment.seq2, cigar)
        print "Gesamtalignment:\n"
        print alignment1.pretty_print_alignment(alignment1.seq1, alignment1.seq2)
        print ""
        tp = alignment1.calculate_intervals(alignment1, verbose)

        if decode:
            tp.rebuild_intervals(tp, tp.seq1, tp.seq2, tp.delta, tp.score, verbose)

    # Ohne CIGAR
    elif args.seq1 and args.seq2 and args.delta and not args.cigar:
        # Nur Trace Points berechnen und speichern
        if decode:
            alignment = alignment1.calculate_alignment(seq1,seq2)
            if verbose:
                print '\nSequenz 1:', seq1
                print 'Sequenz 2:', seq2
                print 'Delta:', delta
                print ""
                print alignment1.pretty_print_alignment(alignment1.seq1, alignment1.seq2)
                print ""
            tp = alignment1.calculate_intervals(alignment,verbose)
        else:
            print '\nSequenz 1:', seq1
            print 'Sequenz 2:', seq2
            print 'Delta:', delta
            print ""
            alignment1 = TracePointAlignment(seq1,seq2,delta) # normales Alignment-Objekt
            alignment = alignment1.calculate_alignment(seq1,seq2) # Alignment mit Gaps und Score
            # Score muss nachgetragen werden, da dieser erst berechnet wird und mit None initialisiert wurde
            alignment1.score = alignment.score

            print "Gesamtalignment:\n"
            print alignment.pretty_print_alignment(alignment.seq1, alignment.seq2)
            print ""

            tp = alignment.calculate_intervals(alignment,verbose) # Alignment mit TracePoints

            print "TP",tp.tp

            if decode:
                tp.rebuild_intervals(tp, tp.seq1, tp.seq2, tp.delta, tp.score, verbose)


    # TODO BAM-File noch nicht vollständig
    elif args.bam:
        bamFile = args.bam
        print "BAM-File:", bamFile
        bamFP = pysam.AlignmentFile(bamFile, "rb")
        # Subprozess in der Konsole
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
            alignment = Alignment.cigar_to_alignment(seq1,seq2,cigar)
            tp = Alignment.calculate_intervals(seq1,seq2,alignment,delta,verbose)
            print "Trace Points:"

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
                                alignment = Alignment.cigar_to_alignment(seq1,seq2,cigar)
                                tp = Alignment.calculate_intervals(seq1,seq2,alignment,delta,verbose)
                                print "Trace Points:"
                        except:
                            print "Fehlerhafter Cigar-String";

    # TODO SAM-File noch nicht vollständig
    elif args.sam:
        samFile = args.sam
        samFP = pysam.Samfile(bamFile, "rb")
        print "SAM-File:", samFile

    # Zufällige Sequenzen
    elif args.random:

        random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
        alignment = TracePointAlignment(seq1,seq2,delta)

        for i in range(0,int(args.random[0])*2,2):
            aln = alignment.calculate_alignment(random_seq_list[i],random_seq_list[i+1])
            if verbose:
                print "\nAlignment zufälliger Sequenzen:"
                print aln.pretty_print_alignment(aln.seq1, aln.seq2)
                print "Score:", aln.score
                print ""
            alignment = TracePointAlignment(random_seq_list[i],random_seq_list[i+1],delta,cigar,aln.score)
            tp = alignment.calculate_intervals(aln,verbose)

            if decode:
                tp.rebuild_intervals(tp, tp.seq1, tp.seq2, tp.delta, tp.score, verbose)
    else:
        sys.stderr.write("Falsche Eingabe der Argumente!")
        sys.exit(1);

if __name__ == "__main__":
    main(sys.argv[1:])

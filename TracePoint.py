# -*- coding: utf-8 -*-

import math
import sys

# Pfad für lokale Submodule
sys.path.append("biopython/")

import string, random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

# Alignment Struktur mit TracePoints
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
    def show_aln(self, seq1, seq2):

        middle = pretty_print= ""
        count = add = 0

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
        pretty_print = "\n" + pos1 + "\n" + seq1 + "\n" + middle + "\n" + seq2 + "\n" + pos2 + "\n"

        return pretty_print

    # berechne Intervalle
    def calculate_intervals(self, alignment, verbose=False):
        count = 0
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

        # Anpassung der Intervalle an Länge der Sequenzen
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
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, letters_in_seq1 - 1, start2, letters_in_seq2 - 1)
                    print self.show_aln(seq1_align[start_seq1:len(seq1_align)], seq2_align[start_seq1:len(seq1_align)])
                seq1_new += str(seq1_align[start_seq1:len(seq1_align)])
                seq2_new += str(seq2_align[start_seq1:len(seq2_align)])

            else:
                # Anzahl der Buchstaben an den Enden der Sequenzen
                rest1 = self.count_indels_letters(seq1_align[count:end_seq1])[0]
                rest2 = self.count_indels_letters(seq2_align[count:end_seq2])[0]

                # Abbruchbedingung
                if rest1 != 0 and rest2 != 0:
                    # first delta letters in seq1_align
                    while self.count_indels_letters(seq1_align[start_seq1:count])[0] != self.delta:
                        count += 1

                indel_count = 0
                indel_count = self.count_indels_letters(seq2_align[start_seq1:count])[1]
                v_id += indel_count
                tp2.append(count - 1 - v_id)
                end = count - 1 - v_id
                indel_count = 0
                if verbose:
                    print "seq1[%d...%d] aligniert mit seq2[%d...%d]" % (start1, k, start2, end)
                    print self.show_aln(seq1_align[start_seq1:count], seq2_align[start_seq1:count])
                    print ""
                seq1_new += str(seq1_align[start_seq1:count])
                seq2_new += str(seq2_align[start_seq1:count])
                start_seq1 = count
                count = start_seq1
                start1 = k + 1
                start2 = start_seq1 - v_id
                end = count - 1 - v_id
        TP_alignment = TracePointAlignment(seq1_new, seq2_new, self.delta, self.score, tp2)

        # Test auf Korrektheit der Sequenzen
        if self.check_alignment(alignment, TP_alignment):
            if verbose:
                print "Konkateniertes Alignment:"
                print self.show_aln(TP_alignment.seq1, TP_alignment.seq2)
            return TP_alignment


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
    def rebuild_intervals(self, tp_alignment, verbose):
        seq1 = tp_alignment.seq1
        seq2 = tp_alignment.seq2
        delta = tp_alignment.delta
        score = tp_alignment.score
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
                    print alignment.show_aln(aln.seq1,aln.seq2)
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
                    print alignment.show_aln(aln.seq1,aln.seq2)
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
                    print alignment.show_aln(aln.seq1, aln.seq2)
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
                    print alignment.show_aln(aln.seq1, aln.seq2)
                    print "Score:",aln.score
                    print ""
                count += 1

        if check_seq1 == seq1 and check_seq2 == seq2:
            if verbose:
                print "Die Konkatenation der Teilalignments ergibt das oben genannte Gesamtalignment!\n"
                print alignment.show_aln(aln1,aln2)
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

    # create TracePoint Alignment
    def create_tp_aln(self, aln, verbose):
        self.calculate_alignment(aln.seq1,aln.seq2)
        if verbose:
            self.show_aln(aln.seq1, aln.seq2)
        tp_aln = self.calculate_intervals(aln,verbose)

        return tp_aln

    # store TracePointAlignment to file
    def store_aln(self, output, mode):

        with open('alignment_compressed.txt', mode) as file_:

            # ';' für einfacheres Splitten
            for item in output:
                file_.write("%s;" % item)
                

    # read TracePointAlignment from file
    def read_and_create_aln(self, aln_file, verbose):
        input_ = aln_seq1 = aln_seq2 = nums = aln_tp = ""
        count = 0

        with open(aln_file, "r") as file_:
            input_ = file_.read()

        # Splitten der einzelnen Elemente
        input_split = input_.split(';')

        # Einlesen der Parameter
        for i in range(0,len(input_split)-1,6):
            aln_seq1 = input_split[i]
            aln_seq2 = input_split[i+1]
            delta = input_split[i+2]
            start_seq1 = input_split[i+3]
            start_seq2 = input_split[i+4]
            aln_tp = input_split[i+5]

            if verbose:
                print "seq1:", aln_seq1
                print "seq2:", aln_seq2
                print "delta:", delta
                print "start_seq1:", start_seq1
                print "start_seq2:", start_seq2
                print "TPs =", aln_tp

            aln = TracePointAlignment(aln_seq1,aln_seq2,delta,start_seq1=start_seq1,start_seq2=start_seq2,tp=tp)
            aln.rebuild_intervals(aln, verbose)



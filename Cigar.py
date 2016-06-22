# -*- coding: utf-8 -*-

import TracePoint

import re

# Alignment Struktur mit Cigar-String
class CigarAlignment(object):

    # Konstruktor
    def __init__(self, seq1, seq2, delta, cigar, start_seq1 = 0, start_seq2 = 0):
        self.seq1 = seq1
        self.seq2 = seq2
        self.delta = delta
        self.cigar = cigar
        self.start_seq1 = start_seq1
        self.start_seq2 = start_seq2

    def pretty_print_cigar_alignment(self, alignment):

        middle = pretty_print = ""
        count = add = 0

        if len(alignment.seq1) < len(alignment.seq2):
            add = len(alignment.seq2) - len(alignment.seq1)
        else:
            add = len(alignment.seq1) - len(alignment.seq2)

        alignment.seq2 += " " * add
        for i in alignment.seq1:
            if i == alignment.seq2[count]:
                middle += '|'
            else:
                middle += ' '
            count += 1

        pos1 = self.print_sequence_positions(alignment.seq1)
        pos2 = self.print_sequence_positions(alignment.seq2)
        pretty_print = "\n" + pos1 + "\n" + alignment.seq1 + "\n" + middle + "\n" + alignment.seq2 + "\n" + pos2 + "\n"
        print pretty_print

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
        alignment = CigarAlignment(seq1_align, seq2_align, self.delta, cigar)

        return alignment

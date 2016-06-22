#!/bin/env python
# -*- coding: utf-8 -*-

import TracePoint
import Cigar

import sys

# Pfad für lokale Submodule
sys.path.append("pysam/")

# für Subprozess in der Konsole
import subprocess
import argparse
import pysam
import string, random
import os


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
                    if r_seq2[position] == 'a':
                        continue
                    r_seq2[position] = 'a'
                elif 0.4 < r_choice < 0.6:
                    # c
                    if r_seq2[position] == 'c':
                        continue
                    r_seq2[position] = 'c'
                elif 0.6 < r_choice < 0.8:
                    # g
                    if r_seq2[position] == 'g':
                        continue
                    r_seq2[position] = 'g'
                elif 0.8 < r_choice <= 1.0:
                    # t
                    if r_seq2[position] == 't':
                        continue
                    r_seq2[position] = 't'
                r_seq2[position] = random.choice(alphabet)
            position += 1
        seq_new1 = "".join(r_seq1)
        seq_new2 = "".join(r_seq2)
        random_seqs.append(seq_new1)
        random_seqs.append(seq_new2)

    return random_seqs


# Berechnung mit Cigar
def with_cigar_calc(cigar_aln, verbose, decode):
    cig_aln = cigar_aln.cigar_to_alignment(cigar_aln.seq1, cigar_aln.seq2, cigar_aln.cigar)
    if verbose:
        print "Cigar-Alignment:"
        cig_aln.pretty_print_cigar_alignment(cig_aln)
    tp_alignment = TracePoint.TracePointAlignment(cig_aln.seq1, cig_aln.seq2, cigar_aln.delta)
    tp = tp_alignment.create_tp_aln(tp_alignment, verbose)


    if decode:
        tp.rebuild_intervals(tp, verbose)


# Berechnung ohne Cigar
def without_cigar_calc(tp_aln, verbose, decode):
    tp = tp_aln.create_tp_aln(tp_aln, verbose)

    if decode:
        tp.rebuild_intervals(tp, verbose)


# Berechnung mit Zufallsgenerator
def random_calc(random_seq_list, delta, verbose, decode):
    for i in range(0, len(random_seq_list), 2):
        tp_alignment = TracePoint.TracePointAlignment(random_seq_list[i], random_seq_list[i + 1], delta)
        if verbose:
            print "\nAlignment zufälliger Sequenzen:"
            print tp_alignment.show_aln(tp_alignment.seq1, tp_alignment.seq2)
            print "Score:", tp_alignment.score
            print ""
        tp = tp_alignment.create_tp_aln(tp_alignment, verbose)

        # Ausgabe in Datei

        output = []
        output.extend((tp.seq1, tp.seq2, delta, tp.start_seq1, tp.start_seq2, (tp.tp)))

        if i == 0:
            # create and write new file
            tp.store_aln(output, 'w')
        elif i == -1:
            # append last output and rebuild if set
            tp.store_aln(output, 'a')
            if decode:
                tp.rebuild_intervals(tp, verbose)

            # TODO noch nicht sinnvoll eingebunden
            # Einlesen aus Datei
            tp.read_and_create_aln("alignment_compressed.txt", verbose)
        else:
            # append other output to existing file
            tp.store_aln(output, 'a')

# Berechnung mit BAM-File
def bam_calc(bamFile, verbose, decode):
    bamFP = pysam.AlignmentFile(bamFile, "rb")
    # Subprozess in der Konsole um mit Samtools die Sequenzen einzulesen
    proc = subprocess.Popen(["samtools view %s| cut -f 10 " % bamFile], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out_split = out.split('\n')
    c = 0

    seq1_list = []
    seq2_list = []

    # Zuweisung von seq1 und seq2
    for s in out_split:
        c += 1
        if c % 2 == 0:
            # seq2 = str(s)
            seq2_list.append(str(s))
        else:
            # seq1 = str(s)
            seq1_list.append(str(s))


    print seq1_list
    print "\n"
    print seq2_list

    # Einlesen von CIGAR String
    cig_type_char = "MIDNSHP"
    cigar = ""

    for read in bamFP:
        if (not (read.is_unmapped)):
            cig_line = read.cigar;
            for (cig_type, cig_length) in cig_line:
                if (cig_type > 6):
                    print "Falsche Cigar-Nummer!"
                    sys.exit(1)
                else:
                    cigar += ("%d" + cig_type_char[cig_type]) % cig_length

    cig_alignment = Cigar.CigarAlignment(seq1, seq2, delta, cigar)
    cig_aln = cig_alignment.cigar_to_alignment(seq1, seq2, cigar)
    if verbose:
        print "Cigar-Alignment:"
        cig_alignment.pretty_print_cigar_alignment(cig_aln)
    tp_alignment = TracePoint.TracePointAlignment(cig_aln.seq1, cig_aln.seq2, delta)
    tp = tp_alignment.create_tp_aln(tp_alignment, verbose)

    if decode:
        tp.rebuild_intervals(tp, verbose)


# Berechnung mit SAM-File
def sam_calc(samFile, verbose, decode):
    return True


def main(argv):
    seq1 = seq2 = cigar = output = ""
    delta = 0

    # Argumente
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq1", "--seq1", help="Die erste Sequenz")
    parser.add_argument("-start_seq1", help="Startposition der ersten Sequenz", type=int, default=0)
    parser.add_argument("-seq2", "--seq2", help="Die zweite Sequenz")
    parser.add_argument("-start_seq2", help="Startposition der zweiten Sequenz", type=int, default=0)
    parser.add_argument("-d", "--delta", help="Delta", type=int)
    group1 = parser.add_mutually_exclusive_group()
    group2 = parser.add_mutually_exclusive_group()
    group1.add_argument("-c", "--cigar", help="Der CIGAR-String")
    group1.add_argument("-b", "--bam", help="Input BAM-File")
    group1.add_argument("-s", "--sam", help="Input SAM-File")
    group1.add_argument("-r", "--random",
                        help="Zufällige Sequenzen generieren mit <Anzahl> <Länge> <Fehlerrate> <Alphabet>",
                        nargs=4)
    group2.add_argument("-e", "--encode", help="Alignment nur mit TracePoints codieren", default=False,
                        action="store_true")
    group2.add_argument("-x", "--decode", help="Aus TracePoints neues Alignment konstruieren", default=False,
                        action="store_true")
    parser.add_argument("-v", "--verbose", help="Ausführlicher Output", default=False, action="store_true")
    args = parser.parse_args()

    seq1 = args.seq1
    start_seq1 = args.start_seq1
    seq2 = args.seq2
    start_seq2 = args.start_seq2
    cigar = args.cigar
    delta = args.delta
    verbose = args.verbose
    decode = args.decode  # initiale Ausgabe des Inputs
    if not args.random:
        print '\nSequenz 1:', seq1
        print 'Sequenz 2:', seq2
        if args.cigar:
            print 'Cigar-String:', cigar
            print 'Delta:', delta
            print ""


    if not args.cigar:
        tp_alignment = TracePoint.TracePointAlignment(seq1, seq2, delta)

    elif args.bam or args.sam:
        pass

    else:
        cigar_alignment = Cigar.CigarAlignment(seq1, seq2, delta, cigar)

    # Mit CIGAR
    if args.seq1 and args.seq2 and args.delta and args.cigar:

        with_cigar_calc(cigar_alignment, verbose, decode)

    # Ohne CIGAR
    elif args.seq1 and args.seq2 and args.delta and not args.cigar:

        without_cigar_calc(tp_alignment, verbose, decode)  

    # Zufällige Sequenzen
    elif args.random:

        random_seq_list = random_sequences(int(args.random[0]), int(args.random[1]), float(args.random[2]), args.random[3])
        random_calc(random_seq_list, delta, verbose, decode)  
    # TODO BAM-File noch nicht vollständig
    elif args.bam:

        bamFile = args.bam
        print "BAM-File:", bamFile
        bam_calc(bamFile, verbose, decode)  

    # TODO SAM-File noch nicht vollständig
    elif args.sam:
        samFile = args.sam
        samFP = pysam.Samfile(bamFile, "rb")
        print "SAM-File:", samFile

    else:
        sys.stderr.write("Falsche Eingabe der Argumente!")
        sys.exit(1);

if __name__ == "__main__":
    main(sys.argv[1:])

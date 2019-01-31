#!/usr/bin/env python

import os

from Bio import SeqIO
from optparse import OptionParser
import Codon_Usage
import gCAI
# import pandas as pd

parser = OptionParser()
parser.add_option("-i", "--dbs", dest="path", default="None",
                  help="[Required] Location of the fasta files to be analyzed")
parser.add_option("-g", "--gc", dest="gc", action="store_true",
                  help="Calculate GC content")
parser.add_option("-f", "--freq", dest="freq", action="store_true",
                  help="Calculate Codon Frequencies")
parser.add_option("-r", "--RCSU", dest="rcsu", action="store_true",
                  help="Calculate RCSU and CAI index")
parser.add_option("-a", "--gCAI", dest="gcai", action="store_true",
                  help="Calculate RCSU and CAI index")
parser.add_option("-c", "--CDS", dest="cds", action="store_true",
                  help="Extract CDS sequences form gb files")


# Parse options into variables
(options, args) = parser.parse_args()

path = options.path
GC = options.gc
FREQ = options.freq
RCSU = options.rcsu
GCAI = options.gcai
CDS_true = options.cds

lsfiles = os.listdir(path)
lsfiles = [i for i in lsfiles if (('fasta' in i) and ("CDS" not in i) and ('cat' not in i))]


def Extract_codons(seq):
    """Extract codon frequency on a given sequence."""
    Codons = {}
    tot = float(len(seq) / 3)
    if len(seq) % 3 == 0 and seq[0:3] == 'ATG':
        for i in range(len(seq) / 3):
            c = seq[i:i + 3].tostring()
            if "N" not in c:
                # print c
                if c not in Codons:
                    Codons[c] = 0
                Codons[c] += 1
    for k in Codons:
        Codons[k] = Codons[k] / tot
    return Codons


def Calculate_GC_Usage(seq):
    """Calculate GC content in sequence (count(G) + count(C) / len(seq) ."""
    G = seq.count('G')
    C = seq.count('C')
    freq = (G + C) / float(len(seq))
    return freq

gcs = {}

for db in lsfiles:
    # for db in ['Drosophila_melanogaster.fasta']:
    print "Doing :", db
    f = open(os.path.join(path, db))
    print "Parsing database ..."
    handle = SeqIO.parse(f, 'gb')
    seqs = []
    for s in handle:
        seqs.append(s)
    print "Recovering CDS ..."
    CDS = []
    for entry in seqs:
        seq = entry.seq
        for f in entry.features:
            if f.type == 'CDS':
                CDS.append(f.extract(seq))
    print "Calculating Codon frequency ..."
    if CDS_true:
        towrite = []
        for i in range(len(CDS)):
            towrite.append(SeqIO.SeqRecord(CDS[i], id=str(i), name=str(i)))
        SeqIO.write(towrite, db.split('.fasta')[0] + '.CDS.fasta', 'fasta')
    if GC or FREQ or RCSU:
        if GC:
            gc = []
        if FREQ:
            codons = []
        for seq in CDS:
            if GC:
                gc.append(Calculate_GC_Usage(seq))
            if FREQ:
                codons.append(Extract_codons(seq))
            if RCSU:
                if not os.path.isfile(os.path.join(path, db.split('.')[0] + '.CDS.fasta')):
                    towrite = []
                    for i in range(len(CDS)):
                        towrite.append(SeqIO.SeqRecord(CDS[i], id=str(i), name=str(i)))
                    SeqIO.write(towrite, os.path.join(path, db.split('.')[0] + '.CDS.fasta'), 'fasta')
                CalcRCSU = Codon_Usage.CodonAdaptationIndex()
                rcsu, cai_index = CalcRCSU.generate_index(os.path.join(path, db.split('.')[0] + '.CDS.fasta'))
                f = open(os.path.join(path, db.split('.')[0] + '.RCSU'), 'w')
                f.write('Codon\tRCSU\n')
                for k in rcsu:
                    f.write('%s\t%s\n' % (k, rcsu[k]))
                f.close()
                f = open(os.path.join(path, db.split('.')[0] + '.CAI_index'), 'w')
                f.write('Codon\tCAI_index\n')
                for k in cai_index:
                    f.write('%s\t%s\n' % (k, cai_index[k]))
                f.close()
    if GCAI:
        if not os.path.isfile(os.path.join(path, db.split('.')[0] + '.CDS.fasta')):
            towrite = []
            for i in range(len(CDS)):
                towrite.append(SeqIO.SeqRecord(CDS[i], id=str(i), name=str(i)))
            SeqIO.write(towrite, os.path.join(path, db.split('.')[0] + '.CDS.fasta'), 'fasta')
        if not os.path.isfile(os.path.join(path, db.split('.')[0] + '.CDS.fasta.highSet')):
            print "Calculating gCAI"
            g = gCAI.gCAI_calculator(os.path.join(path, db.split('.')[0] + '.CDS.fasta'))
            g.Run()
        if not os.path.isfile(os.path.join(path, db.split('.')[0] + '.CodonFreq.tsv')):
            g = gCAI.gCAI_calculator(os.path.join(path, db.split('.')[0] + '.CDS.fasta'))
            g.Calculate_Codon_Frequency()

        print "Saving to dataframe ..."
    # codon_freq = pd.DataFrame(columns=['codon', 'frequency'])

    if GC:
        gcs[db] = gc

    if FREQ:
        codon_freq = open(db.split('.')[0] + ".tsv", 'w')
        codon_freq.write('\t'.join(['codon', 'frequency']) + '\n')
        for c in codons:
            for k in c:
                codon_freq.write('\t'.join([k, str(c[k])]) + '\n')
        codon_freq.close()

if GC:
    gc_freq = open('GC_content' + ".tsv", 'w')
    gc_freq.write('\t'.join(['organism', 'gc_freq']) + '\n')
    for k in gcs:
        for i in gcs[k]:
            gc_freq.write(k.split('.')[0] + '\t' + str(i) + '\n')

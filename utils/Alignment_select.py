#!/usr/bin/env python

# This script is made to rename the HMMER hits turned into a NEWICK tree
# by RAxML back to species, or genus.

import json

from optparse import OptionParser

from Bio import SeqIO

parser = OptionParser()
parser.add_option("-a", "--alignement", dest="alignpath", default="None",
                  help="[Required] Location of the FASTA alignement file")
parser.add_option("-j", "--json", dest="jsonpath", default="None",
                  help="[Required] Location of the JSON file from the HMMER search hits")
parser.add_option("-g", "--goal", dest="goal", default="None",
                  help="[Required] Location of the selected Node names (one per line)")

# Parse options into variables
(options, args) = parser.parse_args()

alignpath = options.alignpath
jsonpath = options.jsonpath
goalpath = options.goal


# load Json File
with open(jsonpath) as data_file:
    metadata = json.load(data_file)

# load Tree
handle = SeqIO.parse(open(alignpath), 'fasta')
Seqs = []
for s in handle:
    Seqs.append(s)

# Parse JSON file into a specie name to accession
name2acc = {}
hits = metadata['results']['hits']
for h in hits:
    if not h['species'] in name2acc:
        name2acc[h['species']] = []
    name2acc[h['species']].append(h['acc'])

goals = []
f = open(goalpath)
for line in f.readlines():
    goals.append(line.strip())

acc2keep = []
for g in goals:
    try:
        for i in name2acc[g]:
            acc2keep.append(i)
    except:
        acc2keep.append(g)

print acc2keep

# grab the sequences which are in the goal file and write into a reduced FASTA file
result = []
for s in Seqs:
    if s.name.split('|')[0] in acc2keep:
        result.append(s)
    else:
        print s.name

print "Number of sequence kept: ", len(result)

f = open(alignpath + '.reduced', 'w')
SeqIO.write(result, f, 'fasta')
f.close()

#!/usr/bin/env python

# This script is made to rename the HMMER hits turned into a NEWICK tree
# by RAxML back to species, or genus.

import json
import sys

from optparse import OptionParser
import progressbar

from Bio import Phylo

parser = OptionParser()
parser.add_option("-t", "--tree", dest="treepath", default="None",
                  help="[Required] Location of the newick tree from RAxML")
parser.add_option("-j", "--taxo", dest="jsonpath", default="None",
                  help="[Required] Location of the taxonomy file from the UNIPROT database")
parser.add_option("-g", "--goal", dest="goal", default="species",
                  help="[Optional] Convert accession numbers to: species, kingdom, description")
parser.add_option("-f", "--format", dest="format", default="nexml",
                  help="[Optional] Output Tree format: newick, nexml, nexus, phyloxml")

# Parse options into variables
(options, args) = parser.parse_args()

treepath = options.treepath
jsonpath = options.jsonpath
Format = options.format
goal = options.goal
if goal not in ['species', 'kingdom', 'description']:
    print "Illegal goal entered"
    sys.exit(1)
else:
    if goal == "kingdom":
        goal = 'kg'
    elif goal == 'description':
        goal = 'desc'

print 'Loading metadata'
metadata = {}
f = open(jsonpath)
lines = f.readlines()
bar = progressbar.ProgressBar()
for i in bar(range(len(lines))):
    line = lines[i]
    s = line.split('\t')
    spec = s[2].split(',')[-1]
    kingdom = s[2].split(',')[0]
    acc = s[0]
    metadata[acc] = kingdom[0] + '|' + spec
    # # load Json File
    # with open(jsonpath) as data_file:
    #     metadata = json.load(data_file)

    # load Tree
trees = Phylo.parse(treepath, 'newick')
Tree = trees.next()

# Modify the name in the tree from Accession numbers to specie name
elements = Tree.get_terminals()
for el in elements:
    name = el.name.split('|')[-1].split('_')[1]
    try:
        specie = metadata[name]
        el.name = specie
    except:
        el.name = name
Phylo.write([Tree], treepath + '.converted', Format)

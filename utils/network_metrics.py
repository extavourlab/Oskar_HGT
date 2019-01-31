#!/usr/bin/env python

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from optparse import OptionParser
from Bio import Phylo

parser = OptionParser()
parser.add_option("-t", "--tree", dest="treepath", default="None",
                  help="[Required] Location of the newick tree from RAxML")
parser.add_option("-n", "--names", dest="namepath", default="None",
                  help="[Required] Location of the file containing nodes of interest for network analysis")


# Parse options into variables
(options, args) = parser.parse_args()

treepath = options.treepath
namepath = options.namepath
if treepath is None or namepath is None:
    print "Invalid options"
    sys.exit(1)

# load Tree
trees = Phylo.parse(treepath, 'newick')
Tree = trees.next()
tree = Phylo.to_networkx(Tree)

# load names
f = open(namepath)
lines = f.readlines()
names = []
for line in lines:
    if line:
        names.append(line.strip())

# Analysis

nodes = tree.nodes()
leaves = []
for node in nodes:
    if node.name is not None:
        leaves.append(node)

paths = {}
paths_lenght = {}

for Ni in range(len(leaves)):
    for Nj in range(len(leaves)):
        if Ni > Nj:
            paths[(Ni, Nj)] = nx.shortest_path(tree, leaves[Ni], leaves[Nj])
            paths_lenght[(Ni, Nj)] = len(paths[(Ni, Nj)]) - 2

names_node = []
for Ni in range(len(leaves)):
    for name in names:
        if name in leaves[Ni].name:
            names_node.append(Ni)

paths_lenght_names = {}
for Ni in names_node:
    for Nj in names_node:
        if Ni > Nj:
            paths_lenght_names[(Ni, Nj)] = paths_lenght[(Ni, Nj)]

names2rest = {}
for k in paths_lenght.keys():
    if k[0] in names_node or k[1] in names_node:
        names2rest[k] = paths_lenght[k]

shortbranch = []
for k in names2rest.keys():
    if names2rest[k] <= 10:
        shortbranch.append(k)

close_species = []
for i in shortbranch:
    if leaves[i[0]].name not in close_species:
        close_species.append(leaves[i[0]].name)
    if leaves[i[1]].name not in close_species:
        close_species.append(leaves[i[1]].name)

print close_species


nodes_dist = np.array([paths_lenght[x] for x in paths_lenght])
nodes_dist_names = np.array([paths_lenght_names[x] for x in paths_lenght_names])
nodes_dist_names2rest = np.array([names2rest[x] for x in names2rest])

sns.distplot(nodes_dist, hist=False, color="g", kde_kws={"shade": True})
sns.distplot(nodes_dist_names, hist=False, color="r", kde_kws={"shade": True})
sns.distplot(nodes_dist_names2rest, hist=False, color="b", kde_kws={"shade": True})
plt.show()

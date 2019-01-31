#!/usr/bin/env python

# import numpy as np
import ete3


class Analysis:
    """Main class for analysis of Tree output."""

    def __init__(self, TreePath, AlignementPath, uniprotTaxonomy):
        """This class takes the path to the Newick Tree, the fasta alignment from which the tree is derived and the path to the parsed uniprot taxonomy."""
        self.TreePath = TreePath
        self.AlignementPath = AlignementPath
        f = open(self.AlignementPath)
        lines = f.readlines()
        out = []
        for line in lines:
            if line[0] == '>':
                out.append(line.split(' ')[0] + '\n')
            else:
                out.append(line)
        f.close()
        f = open(self.AlignementPath, 'w')
        for o in out:
            f.write(o)
        f.close()
        self.tree = ete3.PhyloTree(newick=TreePath, alignment=AlignementPath)
        self.tree.set_species_naming_function(self.parse_sp_name)
        self.uniprot2ncbi = {}
        self.uniprot2species = {}
        self.ncbiID2species = {}
        self.ncbi = ete3.NCBITaxa()
        f = open(uniprotTaxonomy)
        lines = f.readlines()
        for line in lines:
            s = line.strip().split('\t')
            uniprotID = s[0]
            ncbiID = s[1].split(' ')[0]
            specie = s[2].split(',')[-1]
            self.uniprot2ncbi[uniprotID] = ncbiID
            self.uniprot2species[uniprotID] = specie
            self.ncbiID2species[ncbiID] = specie
        self.treeTaxa = []
        leaves = self.tree.get_leaves()
        for leaf in leaves:
            uniprotID = leaf.name.split('|')[0].split('_')[1]
            ncbiID = self.uniprot2ncbi[uniprotID]
            leaf.name = "%s_%s" % (ncbiID, leaf.name.split('|')[0].split('_')[1])
            # leaf.species = sel`f.uniprot2species[uniprotID]
            self.treeTaxa.append(int(ncbiID))
        self.NCBITaxonomy = self.ncbi.get_topology(self.treeTaxa, intermediate_nodes=True)

    def parse_sp_name(self, nodename):
        return self.ncbiID2species[nodename.split("_")[0]]

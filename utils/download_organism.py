#!/usr/bin/env python

# Download all genes for a given specie

from Bio import Entrez
Entrez.email = 'lblondel@g.harvard.edu'


class GeneDownload:

    def __init__(self, organism):
        self.organism = organism
        self.chunk_size = 1000
        self.outfile = open(organism.replace(' ', '_') + '.fasta', 'w')

    def Fetch_Gene_List(self):
        handle = Entrez.esearch(db="nucleotide", term='"%s"[Organism] AND (biomol_mrna[PROP] AND ("201"[SLEN] : "6000"[SLEN]))' % self.organism, retmax=50000)
        record = Entrez.read(handle)
        print "Recovered %s sequences." % record['Count']
        return record['IdList']

    def Download_sequences(self, idList):
        chopped = []
        tmp = []
        for i in range(len(idList)):
            if i % self.chunk_size == 0:
                if i != 0:
                    chopped.append(tmp)
                tmp = []
            tmp.append(idList[i])
        for i in range(len(chopped)):
            print "Downloading %s/%s sequences" % (i * self.chunk_size, len(chopped) * self.chunk_size)
            handle = Entrez.efetch(db="nucleotide", id=",".join(chopped[i]), rettype="gb", retmode="text")
            self.outfile.write(handle.read())

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-i", "--species", dest="species", default="None",
                      help="[Required] Location of the species list file")

    # Parse options into variables
    (options, args) = parser.parse_args()
    species = options.species

    f = open(species)
    lines = f.readlines()
    spec = []
    for l in lines:
        spec.append(l.strip())
    for s in spec:
        print "Doing specie : ", s
        G = GeneDownload(s)
        ids = G.Fetch_Gene_List()
        G.Download_sequences(ids)

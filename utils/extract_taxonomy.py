#!/usr/bin/env python


import sys
from multiprocessing import Pool


def parse(path):
    """Parse uniporot taxonomy file."""
    print "Parsing file: %s" % path
    acc2taxa = {}
    acc2ncbi = {}
    f = open(path)
    line = f.readline()
    tax = []
    while line:
        if line[0:2] == 'ID':
            ID = line.split(' ')[3].split('_')[1]
        if line[0:2] == 'OC':
            [tax.append(i.strip()) for i in line.strip().split('  ')[1].split(';')[:-1]]
        if line[0:2] == 'OX':
            ncbi = line.strip().split('NCBI_TaxID=')[1].split(';')[0]
        if line[0:2] == 'OS':
            name = line.split('  ')[1].strip()
        if line[0:2] == '//':
            # print "Adding %s : %s" % (ID, tax)
            tax.append(name)
            acc2taxa[ID] = tax
            acc2ncbi[ID] = ncbi
            tax = []
        line = f.readline()
    return acc2taxa, acc2ncbi


if __name__ == '__main__':
    P = Pool(8)
    if len(sys.argv) < 2:
        print "Usage: ./extract_taxonomy.py [file1] [file2] [...]"
        sys.exit(1)
    result = P.map(parse, sys.argv[1:])

    print "Writting Results"

    f = open('uniprot_ID_taxa.tsv', 'w')
    for acc2taxa, acc2ncbi in result:
        for k in acc2taxa:
            f.write(k + '\t' + acc2ncbi[k] + '\t' + ','.join(acc2taxa[k]) + '\n')
    f.close()

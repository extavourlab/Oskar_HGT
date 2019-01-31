#!/usr/bin/env python

import json
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-j", "--json", dest="jsonpath", default="None",
                  help="[Required] Location of the JSON file from the HMMER search hits")

# Parse options into variables
(options, args) = parser.parse_args()

jsonpath = options.jsonpath

# load Json File
with open(jsonpath) as data_file:
    metadata = json.load(data_file)

hits = metadata['results']['hits']

sequences = {}
for hit in hits:
    i = 1
    for seq in hit['domains']:
        sequences[hit['acc'] + '|' + str(i) + '|' + seq['cevalue']] = seq['aliaseq'].replace('-', '').upper()

name = hits[0]['domains'][0]['alihmmname']

f = open(name + '.fasta', 'w')

for seq in sequences:
    f.write(">%s\n%s\n" % (seq, sequences[seq]))


from Bio import SeqIO
import numpy as np
from scipy import stats

CodonsDict2 = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
               'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
               'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
               'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
               'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
               'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
               'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
               'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
               'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
               'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
               'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
               'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
               'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

CodonsDict = {'TTT': 0.01, 'TTC': 0.01, 'TTA': 0.01, 'TTG': 0.01, 'CTT': 0.01,
              'CTC': 0.01, 'CTA': 0.01, 'CTG': 0.01, 'ATT': 0.01, 'ATC': 0.01,
              'ATA': 0.01, 'ATG': 0.01, 'GTT': 0.01, 'GTC': 0.01, 'GTA': 0.01,
              'GTG': 0.01, 'TAT': 0.01, 'TAC': 0.01, 'TAA': 0.01, 'TAG': 0.01,
              'CAT': 0.01, 'CAC': 0.01, 'CAA': 0.01, 'CAG': 0.01, 'AAT': 0.01,
              'AAC': 0.01, 'AAA': 0.01, 'AAG': 0.01, 'GAT': 0.01, 'GAC': 0.01,
              'GAA': 0.01, 'GAG': 0.01, 'TCT': 0.01, 'TCC': 0.01, 'TCA': 0.01,
              'TCG': 0.01, 'CCT': 0.01, 'CCC': 0.01, 'CCA': 0.01, 'CCG': 0.01,
              'ACT': 0.01, 'ACC': 0.01, 'ACA': 0.01, 'ACG': 0.01, 'GCT': 0.01,
              'GCC': 0.01, 'GCA': 0.01, 'GCG': 0.01, 'TGT': 0.01, 'TGC': 0.01,
              'TGA': 0.01, 'TGG': 0.01, 'CGT': 0.01, 'CGC': 0.01, 'CGA': 0.01,
              'CGG': 0.01, 'AGT': 0.01, 'AGC': 0.01, 'AGA': 0.01, 'AGG': 0.01,
              'GGT': 0.01, 'GGC': 0.01, 'GGA': 0.01, 'GGG': 0.01}

SynonymousCodons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']
}


class gCAI_calculator:

    def __init__(self, db):
        print "Calculating gCAI for file : ", db
        self.db = db
        self.G = []
        self.S = []
        self.index = {}
        self.G_gCAI = []
        self.S_gCAI = []
        handle = open(db)
        for seq in SeqIO.parse(handle, 'fasta'):
            if len(seq.seq) % 3 == 0:
                skip = False
                alphabet = [i for i in "BDEFHIJKLMNOPQRSUVWXYZ"]
                for let in alphabet:
                    if let in str(seq.seq):
                        skip = True
                if not skip:
                    self.S.append(seq)
        self.total_seq = len(self.S)

    def Load_genome(self, path):
        handle = open(path)
        for seq in SeqIO.parse(handle, 'fasta'):
            if len(seq.seq) % 3 == 0:
                skip = False
                alphabet = [i for i in "BDEFHIJKLMNOPQRSUVWXYZ"]
                for let in alphabet:
                    if let in str(seq.seq):
                        skip = True
                if not skip:
                    self.G.append(seq)

    def Calculate_gCAI(self, seq):
        # stats.gmean()
        gCAI = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            gCAI.append(self.index[codon])
        gCAI = stats.gmean(gCAI)
        # gCAI = 1
        # for i in range(0, len(seq), 3):
        #     codon = seq[i:i + 3]
        #     gCAI *= self.index[codon]
        # gCAI = gCAI**(1.0 / (len(seq) / 3))
        return gCAI

    def Run(self):
        print "Selecting the first Set"
        print "Generating Index"
        self.generate_index()
        print "Calculating gCAI for all sequences"
        self.Calculate_all_gCAI()
        print "Selecting First Set"
        self.set_S()
        print "Starting Iterative search"
        i = 1
        old_S_gCAI = 0
        old_G_gCAI = 0
        self.G_gCAI.append(1)
        while (old_S_gCAI != np.mean(self.S_gCAI)):  # or (max(self.G_gCAI) <= min(self.S_gCAI)):
            print "Step %s" % i
            old_S_gCAI = np.mean(self.S_gCAI)
            old_G_gCAI = np.mean(self.G_gCAI)
            print "Generating Index"
            self.generate_index()
            print "Calculating gCAI for all sequences"
            self.Calculate_all_gCAI()
            print "Selecting Next Set"
            self.set_S()
            print "Old <gCAI> for S = %s" % old_S_gCAI
            print "Old <gCAI> for G = %s" % old_G_gCAI
            print "New <gCAI> for S = %s" % np.mean(self.S_gCAI)
            print "New <gCAI> for G = %s" % np.mean(self.G_gCAI)
            print "max G VS min S = %s VS %s" % (max(self.G_gCAI), min(self.S_gCAI))
            print "Set Lenght = %s (%s)" % (len(self.S), 100 * (len(self.S) / float(len(self.S) + len(self.G))))
            i += 1
            if i > 20:
                break
        SeqIO.write(self.S, self.db + '.highSet', 'fasta')

    def Calculate_Codon_Frequency(self, save_path=None):
        if save_path:
            f = open(save_path, 'w')
        else:
            f = open(self.db.split('.CDS')[0] + '.CodonFreq.tsv', 'w')
            print self.db.split('.CDS')[0] + '.CodonFreq.tsv', 'w'
        col = sorted(CodonsDict.keys())
        for seq in self.G + self.S:
            freq = self._count_codons_freq(seq)
            if freq:
                tmp = []
                for c in col:
                    if c in freq:
                        tmp.append(freq[c])
                    else:
                        tmp.append(0.0)
                # print sum(tmp)
                tmp = [str(i) for i in tmp]
                f.write('\t'.join(tmp) + '\n')
        f.close()

    def Calculate_all_gCAI(self):
        self.G_gCAI = []
        self.S_gCAI = []
        for g in self.G:
            # try:
            self.G_gCAI.append(self.Calculate_gCAI(g.seq))
            # except:
            # print g
        for s in self.S:
            self.S_gCAI.append(self.Calculate_gCAI(s.seq))

    def set_S(self):
        if (len(self.S) / 2.0) > (self.total_seq * 0.01):
            print "Taking 50%"
            toSort = [(self.S[i], self.S_gCAI[i]) for i in range(len(self.S))]
            toSort += [(self.G[i], self.G_gCAI[i]) for i in range(len(self.G))]
            Sorted = sorted(toSort, key=lambda gCAI: gCAI[1], reverse=True)
            print Sorted[0][1], Sorted[100][1], Sorted[-1][1]
            newSet = Sorted[:len(self.S) / 2]
            toG = Sorted[len(self.S) / 2:]
            # if sum([i[1] for i in Sorted[len(Sorted) / 2:]]) > sum([i[1] for i in Sorted[:len(Sorted) / 2]]):
            print min([i[1] for i in newSet]), max([i[1] for i in newSet])
            print min([i[1] for i in toG]), max([i[1] for i in toG])
            newSet = [i[0] for i in newSet]
            toG = [i[0] for i in toG]
            self.S = newSet
            self.G = toG
        else:
            # if len(self.S) - 1 <= int(self.total_seq * 0.01):
            #     toSort = [(self.S[i], self.S_gCAI[i]) for i in range(len(self.S))]
            #     Sorted = sorted(toSort, key=lambda gCAI: gCAI[1], reverse=True)
            #     newSet = [i[0] for i in Sorted[:-1]]
            #     toG = [Sorted[-1][0]]
            #     print "Removing worst gCAI"
            #     print "Number of sequences in S: ", len(newSet)
            #     self.S = newSet
            #     for g in toG:
            #         self.G.append(g)
            # else:
            toSort = [(self.S[i], self.S_gCAI[i]) for i in range(len(self.S))]
            toSort += [(self.G[i], self.G_gCAI[i]) for i in range(len(self.G))]
            Sorted = sorted(toSort, key=lambda gCAI: gCAI[1], reverse=True)
            print "Taking 1%"
            newSet = Sorted[:int(self.total_seq * 0.01)]
            toG = Sorted[int(self.total_seq * 0.01):]
            newSet = [i[0] for i in newSet]
            toG = [i[0] for i in toG]
            self.S = newSet
            self.G = toG

    def Calculate_occurence(self):
        self.occurence = CodonsDict.copy()
        for s in self.S:
            codon_count = self._count_codons_seq(s)
            for codon in codon_count:
                self.occurence[codon] += codon_count[codon]
        for codon in self.occurence:
            self.occurence[codon] = float(self.occurence[codon]) / len(self.S)

    def generate_index(self):
        self.index = CodonsDict.copy()
        # count codon occurrences in the file.
        codon_count = self._count_codons(self.S)
        self.Calculate_occurence()
        # now to calculate the index we first need to sum the number of times
        # synonymous codons were used all together.
        for aa in SynonymousCodons:
            total = 0.0
            rcsu = []  # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
            codons = SynonymousCodons[aa]

            for codon in codons:
                total += codon_count[codon]

            # calculate the RSCU value for each of the codons
            for codon in codons:
                denominator = float(total) / len(codons)
                rcsu.append(codon_count[codon] / denominator)

            # now generate the index W=RCSUi/RCSUmax:
            rcsu_max = max(rcsu)
            for i in range(len(codons)):
                self.index[codons[i]] = self.occurence[codons[i]] * (rcsu[i] / rcsu_max)
                # self.index[codons[i]] = (rcsu[i] / rcsu_max)

    def _count_codons(self, Set):
        # make the codon dictionary local
        codon_count = CodonsDict.copy()
        for seq in Set:
            # Make sure the sequence starts with a start codon and is a multiple of 3
            if seq.seq[0:3] == 'ATG' and len(seq.seq) % 3 == 0:
                # make sure the sequence is lower case
                if str(seq.seq).islower():
                    dna_sequence = str(seq.seq).upper()
                else:
                    dna_sequence = str(seq.seq)
                for i in range(0, len(dna_sequence), 3):
                    codon = dna_sequence[i:i + 3]
                    if codon in codon_count:
                        codon_count[codon] += 1
                    else:
                        print "illegal codon %s in gene: %s" % (codon, seq.id)
        return codon_count

    def _count_codons_seq(self, seq):
        # make the codon dictionary local
        codon_count = CodonsDict2.copy()
        # Make sure the sequence starts with a start codon and is a multiple of 3
        if seq.seq[0:3] == 'ATG' and len(seq.seq) % 3 == 0:
            # make sure the sequence is lower case
            if str(seq.seq).islower():
                dna_sequence = str(seq.seq).upper()
            else:
                dna_sequence = str(seq.seq)
            for i in range(0, len(dna_sequence), 3):
                codon = dna_sequence[i:i + 3]
                if codon in codon_count:
                    codon_count[codon] = 1
                else:
                    print "illegal codon %s in gene: %s" % (codon, seq.id)
        return codon_count

    def _count_codons_freq(self, seq):
        # make the codon dictionary local
        codon_count = CodonsDict2.copy()
        # Make sure the sequence starts with a start codon and is a multiple of 3
        tot = 0.0
        if seq.seq[0:3] == 'ATG' and len(seq.seq) % 3 == 0:
            # make sure the sequence is lower case
            if str(seq.seq).islower():
                dna_sequence = str(seq.seq).upper()
            else:
                dna_sequence = str(seq.seq)
            for i in range(0, len(dna_sequence), 3):
                codon = dna_sequence[i:i + 3]
                if codon in codon_count:
                    tot += 1.0
                    codon_count[codon] += 1
                else:
                    print "illegal codon %s in gene: %s" % (codon, seq.id)
            if tot:
                for c in codon_count:
                    codon_count[c] = codon_count[c] / tot
                return codon_count
            else:
                return None

import random

# DNA Toolbox

Nucleotides = ['A', 'C', 'G', 'T']

RNAcodonDict = {'F': ['UUU', 'UUC'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                'Y': ['UAU', 'UAC'], '*': ['UAA', 'UAG', 'UGA'], 'C': ['UGU', 'UGC'],
                'W': ['UGG'], 'L': ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                'I': ['AUU', 'AUC', 'AUA'],
                'M': ['AUG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'],
                'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                'A': ['GCU', 'GCA', 'GCG', 'GCC'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'],
                'G': ['GGU', 'GGC', 'GGA', 'GGG']}

MonoisotopicMass = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
                    'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
                    'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
                    'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

random_DNA = ''.join([random.choice(Nucleotides) for num in range(30)])


# Check the sequence to make sure it is actually DNA

def ValidateSeq(dna_seq):
    temp_seq = dna_seq.upper()

    for nucl in temp_seq:
        if nucl not in Nucleotides:
            return False
    return temp_seq


# Returns frequency of base pairs
def NucFreq(dna_seq):
    tmpFreq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for nucl in dna_seq:
        tmpFreq[nucl] += 1
    return tmpFreq


# Transcribes from DNA to RNA
def RnaTranscription(dna_seq):
    return dna_seq.replace('T', 'U')


# Complements a given DNA string, same orientation as the given strand
def DnaComplement(dna_seq):
    return dna_seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


# Calculates highest CG percentage in DNA given in a FASTA format
def CGcontent(dna_fasta_format):
    name_dictionary = {}
    percent_dictionary = {}

    splitter = dna_fasta_format.split('>')

    for code in splitter:
        if code != '':
            name_dictionary[code[:13]] = code[13:]

    for key in name_dictionary:
        value = name_dictionary[key]
        CGcount = value.count('C') + value.count('G')
        percent = (CGcount / (len(value))) * 100
        percent_dictionary[key] = percent

    keymax = max(percent_dictionary, key=percent_dictionary.get)

    return keymax, percent_dictionary[keymax]


# Calculate the number of point mutations in a strand(reference other strand)
def HammingDistance(dna_seq, dna_seq2):
    count = 0
    for base in range(len(dna_seq)):
        if dna_seq[base] != dna_seq2[base]:
            count += 1
    return count


# sum([a != b for a, b in zip(s1, s2)])

# Display percent chance of phenotypically dominant offspring given
# k= Homodominant, m=Heterodominant, n=homorecessive, any two can mate
def DomPhenotype(k, m, n):
    population = k + m + n
    total_prob = (k ** 2 + (3 / 4) * (m ** 2) + 2 * m * k + 2 * n * k + m * n - k - m * (3 / 4)) / (
                population * (population - 1))
    return total_prob


# Translation from mRNA to amino acid
def Translation(mRNA):
    codon = ''
    aminoacid = ''
    for base in mRNA:
        if len(codon) < 3:
            codon += base
        if len(codon) == 3:
            for key, value in RNAcodonDict.items():
                if codon in value:
                    aminoacid += key

            codon = ''

    return aminoacid


# Find the positions of the motifs in the DNA
def FindMotif_DNA(dna_seq, dna_motif):
    indices = []
    counter = 0
    for base in dna_seq:
        counter += 1
        if base == dna_motif[0]:
            if dna_seq[counter - 1:counter - 1 + len(dna_motif)] == dna_motif:
                indices.append(counter)

    return indices


# Find the total isotopic mass of a protein
def FindIsotopicMass(amino_acid_seq):
    totalmass = 0
    for aminoacid in amino_acid_seq:
        for key, value in MonoisotopicMass.items():
            if aminoacid == key:
                totalmass += value
    return totalmass


if __name__ == "__main__":
    print(f'DNA Sequence: [{ValidateSeq(random_DNA)}]\n')

    complement = DnaComplement(random_DNA)

    print(f'DNA Complement Strand: [{complement}]\n')

    print(f'Nucleotide Frequency: [{NucFreq(random_DNA)}]\n')

    transcribed = RnaTranscription(random_DNA)

    print(f'RNA Transcription: [{transcribed}]\n')

    translated = Translation(transcribed)

    print(f'mRNA Translation into Amino Acids: [{translated}]\n')

    print(f'Protein Isotopic Mass: [{FindIsotopicMass(translated)}]\n')

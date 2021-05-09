"""
removeRS.py

Remove restriction site from sequence

Modified from code originally written by Will DeLoache.
"""

from Bio.Seq import Seq
from Bio.Restriction import Restriction, FormattedSeq

GeneticCode = {'CTT': 'L', 'TAG': '*', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
               'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
               'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
               'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
               'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P',
               'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R',
               'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F',
               'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TAA': '*',
               'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L',
               'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A',
               'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V',
               'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L',
               'CGT': 'R', 'TGG': 'W', 'ATG': 'M', 'CGC': 'R'}

HsCodonUsage = {'A': [('GCT', 0.27), ('GCA', 0.23), ('GCC', 0.40), ('GCG', 0.11)],
                'C': [('TGT', 0.46), ('TGC', 0.54)],
                'E': [('GAA', 0.42), ('GAG', 0.58)],
                'D': [('GAT', 0.46), ('GAC', 0.54)],
                'G': [('GGT', 0.16), ('GGA', 0.25), ('GGC', 0.34), ('GGG', 0.25)],
                'F': [('TTT', 0.46), ('TTC', 0.54)],
                'I': [('ATT', 0.36), ('ATA', 0.17), ('ATC', 0.47)],
                'H': [('CAT', 0.42), ('CAC', 0.58)],
                'K': [('AAA', 0.43), ('AAG', 0.57)],
                '*': [('TAA', 0.30), ('TGA', 0.47), ('TAG', 0.24)],
                'M': [('ATG', 1.0)],
                'L': [('TTG', 0.13), ('TTA', 0.08), ('CTA', 0.07), ('CTT', 0.13), ('CTG', 0.40), ('CTC', 0.20)],
                'N': [('AAT', 0.47), ('AAC', 0.53)],
                'Q': [('CAA', 0.27), ('CAG', 0.73)],
                'P': [('CCA', 0.28), ('CCT', 0.29), ('CCC', 0.32), ('CCG', 0.11)],
                'S': [('TCT', 0.19), ('TCA', 0.15), ('AGT', 0.15), ('TCC', 0.22), ('AGC', 0.24), ('TCG', 0.05)],
                'R': [('AGA', 0.21), ('AGG', 0.21), ('CGT', 0.08), ('CGA', 0.11), ('CGC', 0.18), ('CGG', 0.20)],
                'T': [('ACT', 0.25), ('ACA', 0.28), ('ACC', 0.36), ('ACG', 0.11)],
                'W': [('TGG', 1.0)],
                'V': [('GTT', 0.18), ('GTA', 0.12), ('GTC', 0.24), ('GTG', 0.46)],
                'Y': [('TAT', 0.44), ('TAC', 0.56)]}


def findSite(seq, enzyme):
    site = enzyme.site
    if site.find("N") != -1:
        raise Exception("Degenerate enzymes are not supported")
    u_seq = seq.upper()
    loc1 = u_seq.find(site)
    loc2 = u_seq.find(_reversecomplement(site))
    if loc1 == -1 and loc2 == -1:
        return False
    if loc1 == -1:
        return loc2
    if loc2 == -1:
        return loc1
    else:
        return min(loc1, loc2)


# returns false if the codons aren't exactly 1 bp different
# returns index of mismatched base if only 2 matches
def singleBPmutation(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    matches = 0
    mismatchedBase = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matches += 1
        else:
            mismatchedBase = i
    if matches == 2:
        return mismatchedBase
    else:
        return False


def silentMutate(seq, leftIndex, rightIndex, enzyme_list=[]):
    firstCodonIndex = leftIndex - (leftIndex % 3)
    numCodons = int(rightIndex / 3) - int(leftIndex / 3) + rightIndex % 3
    possCodons = []
    # i holds the codon number
    for i in range(numCodons):
        currCodon = seq[firstCodonIndex + i * 3:firstCodonIndex + i * 3 + 3]
        if len(currCodon) == 3:
            allPossCodons = HsCodonUsage[GeneticCode[currCodon.upper()]]
            for codon in allPossCodons:
                # p holds the position within the codon of the mutated base (should always be 2 as written currently)
                p = singleBPmutation(codon[0], currCodon)
                if p and leftIndex <= p + i * 3 + firstCodonIndex < rightIndex:
                    possCodons.append([codon[0], codon[1], i, p])
    possCodons.sort(key=lambda x: x[1])
    if len(possCodons) < 1:
        raise Exception("Couldn't find a base to mutate silently.")
    successfullyMutated = False
    while not successfullyMutated:
        newCodon = possCodons.pop()

        oldBase = seq[
                  firstCodonIndex + newCodon[2] * 3 + newCodon[3]:firstCodonIndex + newCodon[2] * 3 + newCodon[3] + 1]
        newBase = ""
        if oldBase.islower():
            newBase = newCodon[0][newCodon[3]].upper()
        else:
            newBase = newCodon[0][newCodon[3]].lower()

        mutationIndex = firstCodonIndex + newCodon[2] * 3 + newCodon[3]
        leftBase = seq[mutationIndex - 1:mutationIndex]
        rightBase = seq[mutationIndex + 1:mutationIndex + 2]
        front = seq[:mutationIndex]
        if (leftBase.islower() and newBase.islower()) or (leftBase.isupper() and newBase.isupper()):
            front = front.swapcase()
        back = seq[mutationIndex + 1:]
        if (rightBase.islower() and newBase.islower()) or (rightBase.isupper() and newBase.isupper()):
            back = back.swapcase()
        newSeq = front + newBase + back

        introducedNewSite = False
        for enzyme_name in enzyme_list:
            enzyme = getattr(Restriction, enzyme_name)
            orig = FormattedSeq(Seq(seq))
            new = FormattedSeq(Seq(newSeq))
            if len(enzyme.search(new)) > len(enzyme.search(orig)):
                introducedNewSite = True
        if not introducedNewSite or len(possCodons) < 1:
            successfullyMutated = True

    s1 = Seq(seq)
    s2 = Seq(newSeq)
    if str(s1.translate()) != str(s2.translate()):
        raise Exception("Error: The attempted silent mutation wasn't silent!!!")
    return newSeq


def removeRS(seq, enzymes):
    removeRS_tuples = list()
    mutated = True
    while mutated:
        mutated = False
        for enzyme_name in enzymes:
            enzyme = getattr(Restriction, enzyme_name)
            while findSite(seq, enzyme):
                leftIndex = findSite(seq.lower(), enzyme)
                rightIndex = leftIndex + len(enzyme.site)
                seq = silentMutate(seq, leftIndex, rightIndex, enzyme_list=enzymes)
                removeRS_tuples.append((leftIndex, rightIndex))
                mutated = True
    return seq, removeRS_tuples


def _reversecomplement(sequence):
    """Return the reverse complement of the dna string."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    reverse_complement_sequence = ""
    sequence_list = list(sequence)
    sequence_list.reverse()
    for letter in sequence_list:
        reverse_complement_sequence += complement[letter.upper()]
    return reverse_complement_sequence

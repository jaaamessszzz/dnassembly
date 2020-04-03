#! /usr/bin/env python3
from itertools import tee

dna_basepairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Thank you itertools cookbook!
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def cycle_in_frames(iterable, frame=3):
    """
    Cycle through iterable and yield all possible frames of size frame. Treats iterable as a circular sequence.
    Janky but it works...
    :return:
    """
    assert frame < len(iterable)
    cycle_iterable = iterable * 2

    for cycle_count, derp in enumerate(cycle_iterable):
        if cycle_count < len(iterable):
            yield cycle_iterable[cycle_count: cycle_count + frame]
            cycle_count += 1


def reverse_complement(seq):
    "Get reverse complement of DNA sequence"
    return ''.join(reversed([dna_basepairs[a] if a in dna_basepairs.keys() else a for a in seq]))


# --- Manually implemented b/c Biopython is a mess --- #

# E. coli codon table: http://www.sci.sdsu.edu/~smaloy/MicrobialGenetics/topics/in-vitro-genetics/codon-usage.html
# Listed in order of frequency
res_to_codons = {'A': ['GCG', 'GCC', 'GCA', 'GCT'],
                 'C': ['TGC', 'TGT'],
                 'D': ['GAT', 'GAC'],
                 'E': ['GAA', 'GAG'],
                 'F': ['TTT', 'TTC'],
                 'G': ['GGC', 'GGT', 'GGG', 'GGA'],
                 'H': ['CAT', 'CAC'],
                 'I': ['ATT', 'ATC', 'ATA'],
                 'K': ['AAA', 'AAG'],
                 'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
                 'M': ['ATG'],
                 'N': ['AAC', 'AAT'],
                 'P': ['CCG', 'CCA', 'CCT', 'CCC'],
                 'Q': ['CAG', 'CAA'],
                 'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
                 'S': ['AGC', 'TCT', 'TCC', 'TCG', 'AGT', 'TCA'],
                 'T': ['ACC', 'ACA', 'ACG', 'ACT'],
                 'V': ['GTG', 'GTT', 'GTC', 'GTA'],
                 'W': ['TGG'],
                 'Y': ['TAT', 'TAC'],
                 '*': ['TAA', 'TGA', 'TAG'],
                 }

codon_to_res = {codon: res for res, codons in res_to_codons.items() for codon in codons}

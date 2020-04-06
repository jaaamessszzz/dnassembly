#! /usr/bin/env python3

"""
Generate primers from an input DNA sequence
"""
import re

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

from ..utils import reverse_complement

def create_amplifiction_primers(sequence, prefix='', suffix='', target_tm=60, primer=50, Na=50, K=None, Mg=None, dNTPs=None, Tris=None):
    """
    Create amplification primers for a given sequence

    :param sequence: DNA object or string
    :param prefix:
    :param suffix:
    :param moclo_standardize: DNA object or string
    :return: A list of DNA primers
    """

    # Make sure sequence is DNA (ATCG)
    if re.fullmatch('[ATCG]+', sequence.upper()) is None:
        raise Exception('Input sequence is not valid DNA!')
    if prefix != '' and re.fullmatch('[ATCG]+', prefix.upper()) is None:
        raise Exception('Input prefix is not valid DNA!')
    if suffix != '' and re.fullmatch('[ATCG]+', suffix.upper()) is None:
        raise Exception('Input suffix is not valid DNA!')

    # Forward
    forward_primer = generate_primer(sequence, prefix=prefix, target_tm=target_tm, primer=primer, Na=Na, K=K, Mg=Mg, dNTPs=dNTPs, Tris=Tris)

    # Reverse
    sequence_revcomp = reverse_complement(sequence)
    suffix_revcomp = reverse_complement(suffix)
    reverse_primer = generate_primer(sequence_revcomp, prefix=suffix_revcomp, target_tm=target_tm, primer=primer, Na=Na, K=K, Mg=Mg, dNTPs=dNTPs, Tris=Tris)

    return forward_primer, reverse_primer


def generate_primer(sequence, prefix='', target_tm=60, primer=50, Na=50, K=None, Mg=None, dNTPs=None, Tris=None):
    """
    Create an amplification primer for a given sequence

    :param sequence: DNA string
    :return: A list of DNA primers
    """

    # Make sure sequence is DNA (ATCG)
    if re.fullmatch('[ATCG]+', sequence.upper()) is None:
        raise Exception('Input sequence is not valid DNA!')
    if prefix != '' and re.fullmatch('[ATCG]+', prefix.upper()) is None:
        raise Exception('Input prefix is not valid DNA!')

    target_primer = None

    for primer_length in range(len(sequence) + 1):
        primer_TM = mt.Tm_NN(Seq(sequence[:primer_length]),
                             nn_table=mt.DNA_NN2,
                             dnac1=primer / 2,  # nM Primers / 2
                             dnac2=primer / 2,  # nM Primers / 2
                             selfcomp=False,
                             Na=Na,  # mM
                             K=K or 0,  # mM
                             Tris=Tris or 0,  # mM
                             Mg=Mg or 0,  # mM
                             dNTPs=dNTPs or 0,
                             saltcorr=5)
        if primer_TM >= target_tm:
            target_primer = sequence[:primer_length]
            break

    if target_primer is None:
        print('A primer with the desired TM could not be generated for the input sequence!')
        return None

    return prefix + target_primer


def create_ligation_primers(sequence):
    """
    Generate primers for ligation block(s) in Golden Gate Assembly, assuming the input sequence includes overhangs
    If input sequence is >60bp, primers for multiple ligation blocks will be returned
    Reference: https://www.idtdna.com/pages/education/decoded/article/annealing-oligonucleotides

    :param sequence:
    :return:
    """

    if len(sequence) > 60:
        # Break into multiple segments made to be as close to 60bp as possible
        pass
    else:
        forward_primer = sequence[:-4]
        reverse_primer = reverse_complement(sequence[4:])



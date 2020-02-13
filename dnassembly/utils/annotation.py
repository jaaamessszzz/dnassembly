#! /usr/bin/env python3

import re

from Bio.Restriction import BsaI, BsmBI

from ..reactions.moclo import PartOrder, CassetteOrder
from .utils import reverse_complement

bsai_annotation_f = {'CCCT': '1',
                     'AACG': '2',
                     'TATG': '3a',
                     'TTCT': '3b',
                     'ATCC': '4a',
                     'TGGC': '4b',
                     'GCTG': '5',
                     'TACA': '6',
                     'GAGT': '7',
                     'CCGA': '8a',
                     'CAAT': '8b'
                     }

bsai_annotation_r = {'AACG': '1',
                     'TATG': '2',
                     'TTCT': '3a',
                     'ATCC': '3b',
                     'TGGC': '4a',
                     'GCTG': '4b',
                     'TACA': '5',
                     'GAGT': '6',
                     'CCGA': '7',
                     'CAAT': '8a',
                     'CCCT': '8b'
                     }

bsmbi_annotation_f = {'CTGA': 'LS-1',
                      'CCAA': '1-2',
                      'GATG': '2-3',
                      'GTTC': '3-4',
                      'GGTA': '4-5',
                      'AAGT': '5-RE',
                      'AGCA': 'RE-LS'}

bsmbi_annotation_r = {'CCAA': 'LS-1',
                      'GATG': '1-2',
                      'GTTC': '2-3',
                      'GGTA': '3-4',
                      'AAGT': '4-5',
                      'AGCA': '5-RE',
                      'CTGA': 'RE-LS'}


def annotate_moclo(dna, annotate='part'):
    """
    Read the dna sequence and report any parts contained within. Parts are flanked by BsaI sequences.

    Returns None if no part is found, otherwise a list of parts in order of 5'->3'
    """
    
    if annotate == 'part':
        OrderLinkedList = PartOrder()
        RxnEnzyme = BsaI
        annotation_f = bsai_annotation_f
        annotation_r = bsai_annotation_r
    elif annotate == 'cassette':
        OrderLinkedList = CassetteOrder()
        RxnEnzyme = BsmBI
        annotation_f = bsmbi_annotation_f
        annotation_r = bsmbi_annotation_r
    else:
        return Exception("Annotation type must be 'part' or 'cassette'!")

    rxnsite_f = re.finditer(f'{RxnEnzyme.site}', dna)
    rxnsite_r = re.finditer(f'{reverse_complement(RxnEnzyme.site)}', dna)

    match_f = 0
    match_r = 0

    # Get forward sticky end
    for match in rxnsite_f:
        stick_end_f = dna[match.start() + RxnEnzyme.fst5: match.end() + RxnEnzyme.fst3]
        match_f += 1

    # Get reverse sticky end
    for match in rxnsite_r:
        stick_end_r = dna[match.start() - RxnEnzyme.fst3: match.end() - RxnEnzyme.fst5]
        match_r += 1

    if match_r != 1 or match_f != 1:
        return False

    current_part = OrderLinkedList.head
    print(current_part.part)

    part_end = False
    add_parts = False
    part_list = []

    while part_end is False:
        forward_part = annotation_f[stick_end_f]
        reverse_part = annotation_r[stick_end_r]

        if current_part.part == forward_part:
            add_parts = True

        print(add_parts)
        print(current_part.part, forward_part)
        print(current_part.part, reverse_part)

        if current_part.part == reverse_part and add_parts:
            part_end = True

        if add_parts:
            part_list.append(current_part.part)

        current_part = current_part.next
        print(current_part.part)

    return part_list


# def annotate_moclo_parts(dna):
#     """
#     Read the dna sequence and report any parts contained within. Parts are flanked by BsaI sequences.
# 
#     Returns None if no part is found, otherwise a list of parts in order of 5'->3'
#     """
# 
#     PartOrderLinkedList = PartOrder()
# 
#     rxnsite_f = re.finditer(f'{BsaI.site}', dna)
#     rxnsite_r = re.finditer(f'{reverse_complement(BsaI.site)}', dna)
# 
#     match_f = 0
#     match_r = 0
# 
#     # Get forward sticky end
#     for match in rxnsite_f:
#         stick_end_f = dna[match.start() + BsaI.fst5: match.end() + BsaI.fst3]
#         match_f += 1
# 
#     # Get reverse sticky end
#     for match in rxnsite_r:
#         stick_end_r = dna[match.start() - BsaI.fst3: match.end() - BsaI.fst5]
#         match_r += 1
# 
#     if match_r != 1 or match_f != 1:
#         return False
# 
#     current_part = PartOrderLinkedList.head
#     print(current_part.part)
# 
#     part_end = False
#     add_parts = False
#     part_list = []
# 
#     while part_end is False:
#         forward_part = bsai_annotation_f[stick_end_f]
#         reverse_part = bsai_annotation_r[stick_end_r]
# 
#         if current_part.part == forward_part:
#             add_parts = True
# 
#         print(add_parts)
#         print(current_part.part, forward_part)
#         print(current_part.part, reverse_part)
# 
#         if current_part.part == reverse_part and add_parts:
#             part_end = True
# 
#         if add_parts:
#             part_list.append(current_part.part)
# 
#         current_part = current_part.next
#         print(current_part.part)
# 
#     return part_list

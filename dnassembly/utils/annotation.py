#! /usr/bin/env python3

import re

from Bio.Restriction import BsaI, BsmBI

from ..reactions.moclo import PartOrder, CassetteOrder
from .utils import reverse_complement

def annotate_moclo(dna, annotate='part'):
    """
    Read the dna sequence and report any parts contained within. Parts are flanked by BsaI sequences.

    Returns None if no part is found, otherwise a list of parts in order of 5'->3'
    """
    
    if annotate == 'part':
        OrderLinkedList = PartOrder()
        RxnEnzyme = BsaI
        annotation_f = OrderLinkedList.bsai_annotation_f
        annotation_r = OrderLinkedList.bsai_annotation_r
    elif annotate == 'cassette':
        OrderLinkedList = CassetteOrder()
        RxnEnzyme = BsmBI
        annotation_f = OrderLinkedList.bsmbi_annotation_f
        annotation_r = OrderLinkedList.bsmbi_annotation_r
    else:
        return Exception("Annotation type must be 'part' or 'cassette'!")

    rxnsite_f = re.finditer(f'{RxnEnzyme.site}', dna.upper())
    rxnsite_r = re.finditer(f'{reverse_complement(RxnEnzyme.site)}', dna.upper())

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

    forward_part = annotation_f.get(stick_end_f.upper())
    reverse_part = annotation_r.get(stick_end_r.upper())

    if forward_part is None or reverse_part is None:
        return False

    current_part = OrderLinkedList.head

    part_end = False
    add_parts = False
    part_list = []

    while part_end is False:

        if current_part.part == forward_part:
            add_parts = True

        if current_part.part == reverse_part and add_parts:
            part_end = True

        if add_parts:
            part_list.append(current_part.part)

        current_part = current_part.next

    if len(part_list) == len(OrderLinkedList.parts):
        return False

    return part_list

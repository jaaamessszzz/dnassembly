#! /usr/bin/env python3

import re
from Bio.Restriction import BbsI, BsmBI
from typing import List
from enum import Enum
from .utils import reverse_complement
from ..reactions.moclo import PartOrder, CassetteOrder

class AnnotationType(Enum):
    PART = 'part'
    CASSETTE = 'cassette'

def annotate_moclo(dna: str, annotate: str = 'part') -> List[str]:
    """Read the dna sequence and report any parts contained within.
    Parts are flanked by BsaI sequences.

    Args:
        dna (str): dna sequence
        annotate (str, optional): Annotation type. Defaults to 'part'.
    
    TODO - revise this return documentation
    Returns:
        List[str]: Returns None if no part is found, otherwise a list of parts in order of 5'->3'
    """

    #TODO - stick bbsi_annotation_f and bbsi_annotation_r in config file?

    if annotate == 'part':
        OrderLinkedList = PartOrder()
        RxnEnzyme = BbsI
        annotation_f = OrderLinkedList.bbsi_annotation_f
        annotation_r = OrderLinkedList.bbsi_annotation_r
    elif annotate == 'cassette':
        OrderLinkedList = CassetteOrder()
        RxnEnzyme = BsmBI
        annotation_f = OrderLinkedList.bsmbi_annotation_f
        annotation_r = OrderLinkedList.bsmbi_annotation_r
    else:
        raise Exception("Annotation type must be 'part' or 'cassette'!")

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

    #TODO - should convert this return type to something else
    if len(part_list) == len(OrderLinkedList.parts):
        return False

    return part_list

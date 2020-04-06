#! /usr/bin/env python3

import re
from enum import Enum

from Bio.Restriction import BsaI, BsmBI

from ..dna import Part
from .goldengate import GoldenGate
from ..utils.utils import res_to_codons, codon_to_res


# --- Circular double linked lists for Part/Cassette Annotation --- #

class PartNode(object):
    """
    Node in circular double linked list
    """
    def __init__(self, part):
        self.part = part
        self.previous = None
        self.next = None


class CircularOrder(object):
    """
    Base class for circular double linked list
    """
    def __init__(self):
        self.head = None
        self.size = 0
        self.assemblelinkedlist()

    def add(self, part):
        new_part = PartNode(part)

        if self.head is None:
            self.head = new_part
            self.head.previous = self.head
            self.head.next = self.head

        else:
            tail = self.head.previous
            tail.next = new_part
            new_part.previous = tail
            new_part.next = self.head
            self.head.previous = new_part

        self.size += 1

    def assemblelinkedlist(self):
        for part in self.parts:
            self.add(part)


class PartOrder(CircularOrder):
    """
    Part plasmid annotation
    """

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

    parts = [part for part in bsai_annotation_f.values()]

    def __init__(self):
        super(PartOrder, self).__init__()


class CassetteOrder(CircularOrder):
    """
    Cassette plasmid annotation
    """

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

    parts = [part for part in bsmbi_annotation_f.values()]

    def __init__(self):
        super(CassetteOrder, self).__init__()


class MoCloAssemblyType(Enum):
    PART = BsmBI
    CASSETTE = BsaI
    MULTICASSETTE = BsmBI


# --- MoClo Cloning Function --- #

class ModularCloning(GoldenGate):
    """
    The MoClo system is just Golden Gate assembly with extra rules
    """

    def __init__(self, parts, type):
        super(ModularCloning, self).__init__(parts)
        self.restriction_enzyme_list = type
        self.assembly_type = type

    @property
    def restriction_enzyme_list(self):
        return self._restriction_enzyme_list

    @restriction_enzyme_list.setter
    def restriction_enzyme_list(self, type):
        if isinstance(type, MoCloAssemblyType):
            self._restriction_enzyme_list = [type.value]
        else:
            raise AssemblyTypeException('Assembly type must be set using the MoCloAssemblyType class!')


# --- MoClo Related Functions --- #

def MoCloPartFromSequence(sequence, part_5, part_3, description=None, standardize=True, create_instructions=False, remove_bsai=True, remove_bsmbi=True, remove_noti=True):
    """
    Create a MoClo compatible part from an arbitrary sequence
    Checks for BsaI/BsmBI/NotI restriction sites
    Verifies Part 3 definitions are in frame

    Some observations:
    * BsaI/BsmBI/NotI sites can only share max 2bp overlap (BsmBI/NotI)

    :param sequence: string for a DNA sequence, this will be immediately converted to a DNA object
    :param part_5: Left overhang part definition, must be a part defined in PartOrder
    :param part_3: Right overhang part definition, must be a part defined in PartOrder
    :param description: part description
    :returns Part object for the input sequence with the specified overhangs
    """
    if part_5 not in PartOrder.parts or part_3 not in PartOrder.parts:
        raise Exception('Invalid part definitions were passed!')

    # Make sure Part 3 definitions are in frame
    if '3' in part_5 and '3' in part_3:
        if len(sequence) % 3 != 0:
            raise Exception('Part 3 coding sequence definitions must be in frame!')

    # "GG" is appended to any Part 3 sequence to abide by the GS linker definition in the YTK
    standardize_5 = {'3a': 'GG',
                     '3b': 'GG',
                     }

    standardize_3 = {'4a': 'TAA',
                     }

    # Add part overhangs
    overhang_5 = {v:k for k,v in PartOrder.bsai_annotation_f.items()}
    overhang_3 = {v:k for k,v in PartOrder.bsai_annotation_r.items()}

    standardized_5 = standardize_5[part_5] if (part_5 in standardize_5.keys() and standardize) else ""
    standardized_3 = standardize_3[part_3] if (part_3 in standardize_3.keys() and standardize) else ""

    prefix = f'GCATCGTCTCATCGGTCTCA{overhang_5[part_5]}{standardized_5}'
    suffix = f'{standardized_3}{overhang_3[part_3]}TGAGACCTGAGACGGCAT'

    # Create regex pattern for restriction sites
    rxn_regex = []
    if remove_bsai:
        rxn_regex += ['GGTCTC', 'GAGACC']
    if remove_bsmbi:
        rxn_regex += ['CGTCTC', 'GAGACG']
    if remove_noti:
        rxn_regex += ['GCGGCCGC', 'CGCCGGCG']

    regex_pattern = f"?=({'|'.join(rxn_regex)})"

    if create_instructions:
        create_assembly_instructions(sequence, prefix, suffix, part_5, part_3, remove_bsai=remove_bsai, remove_bsmbi=remove_bsmbi, remove_noti=remove_noti)
    else:
        final_sequence = prefix + sequence + suffix
        sequence_DNA = Part(final_sequence, description=description)
        if len(re.findall(regex_pattern, sequence)) > 0:
            raise Exception('There are BsaI/BsmBI/NotI restriction sites it your part definition! Please remove them.')
        return sequence_DNA


# todo: figure this out...
def create_assembly_instructions(sequence, part_5, part_3, prefix='', suffix='', remove_bsai=True, remove_bsmbi=True, remove_noti=True):
    """Create assembly instructions for a part"""
    # Create regex pattern for restriction sites
    rxn_regex = []
    if remove_bsai:
        rxn_regex += ['GGTCTC', 'GAGACC']
    if remove_bsmbi:
        rxn_regex += ['CGTCTC', 'GAGACG']
    if remove_noti:
        rxn_regex += ['GCGGCCGC', 'CGCCGGCG']

    regex_pattern = f"?=({'|'.join(rxn_regex)})"

    # Identify potential restriction sites in left_arm+sequence+right_arm
    matches = re.finditer(regex_pattern, sequence)
    match_results = [match.start() for match in matches]

    # separate sequence into codons
    sequence_codons = [sequence[codon:codon + 3].upper() for codon in range(0, len(sequence), 3)]

    # Track codons
    codon_substitutions = {}

    def make_codon_substitution(sequence_index, bad_codons=None):
        """Substitute a codon in sequence_codons"""
        codon_index = sequence_index // 3
        left_codon_slice = codon_index - 2 if codon_index - 2 >= 0 else 0
        right_codon_slice = codon_index + 2 if codon_index + 2 < len(sequence_codons) else len(sequence_codons)
        local_sequence_codons = sequence_codons[left_codon_slice:right_codon_slice]

        # Check if a previous substitution fixed this already, otherwise fix
        if len(re.findall(regex_pattern, ''.join(local_sequence_codons))) > 0:

            substitution_successful = False

            for index_step in range(3):
                current_codon_index = codon_index + index_step
                codon_to_swap = sequence_codons[current_codon_index]

                codons_to_exclude = [codon_to_swap]
                if bad_codons is not None and type(bad_codons) == list:
                    codons_to_exclude += bad_codons
                codons_to_try = [codon for codon in codon_to_res[res_to_codons[codon_to_swap]] if
                                 codon not in codons_to_exclude]

                # Swap out a single codon per match
                for codon in codons_to_try:
                    test_local_sequence = local_sequence_codons
                    test_local_sequence[2 + index_step] = codon
                    local_seqeunce = ''.join(test_local_sequence)
                    local_matches = re.findall(regex_pattern, local_seqeunce)

                    if len(local_matches) == 0:
                        codon_substitutions[current_codon_index] = {}
                        codon_substitutions[current_codon_index]['original'] = sequence_codons[current_codon_index]
                        codon_substitutions[current_codon_index]['new'] = codon

                        # Make substitution
                        sequence_codons[current_codon_index] = codon
                        substitution_successful = True
                        break

                if substitution_successful:
                    break

            if not substitution_successful:
                raise Exception("A substitution wasn't made somehow...")

    # Identify codons that need to be changed
    # Max of two consecutive codons needed to guarantee removal of rxn site
    for match in match_results:
        make_codon_substitution(match)

    # todo: check for potential restriction sites formed when left/right arms are added to seqeunce

    # Pick overhang sequences
    overhang_5 = {v:k for k,v in PartOrder.bsai_annotation_f.items()}
    overhang_3 = {v:k for k,v in PartOrder.bsai_annotation_r.items()}

    codon_indicies = sorted(list(codon_substitutions.keys()))
    used_overhang_sequences = [overhang_5[part_5], overhang_3[part_3]]
    overhang_start_indicies = []

    def try_overhangs(sequence_index):
        """Choose one of two possible overhangs"""
        if sequence[sequence_index:sequence_index + 4] not in used_overhang_sequences:
            used_overhang_sequences.append(sequence[sequence_index:sequence_index + 4])
            overhang_start_indicies.append(sequence_index)
        elif sequence[sequence_index - 1:sequence_index + 3] not in used_overhang_sequences:
            used_overhang_sequences.append(sequence[sequence_index - 1:sequence_index + 3])
            overhang_start_indicies.append(sequence_index - 1)
        else:
            make_codon_substitution(sequence_index, bad_codons=[codon_substitutions[sequence_index // 3]['original'],
                                                                codon_substitutions[sequence_index // 3]['new']])
            try_overhangs(sequence_index)  # This is sketch...

    # Assign codon overhangs, bad overhangs will be ignored in block generation step
    for count, codon_index in enumerate(codon_indicies):
        # Edge case: first codon in sequence, ligation block if next substitution is close else gBlock
        if codon_index == 0:
            overhang_start_indicies.append(0)
        # Edge case: last codon in sequence, ligation block if previous substitution is close else gBlock
        if codon_index == len(sequence_codons):
            overhang_start_indicies.append(len(sequence - 4))

        sequence_index = codon_index * 3
        try_overhangs(sequence_index)

    """
    If len(codon_substitutions) > 3, just order a gBlock
    Default to PCR
    If substitutions were made <100 nt apart, make ligation blocks
    """

    # todo: get new potentially altered sequence


class AssemblyTypeException(Exception):
    pass

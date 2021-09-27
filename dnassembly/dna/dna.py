#! /usr/bin/env python3

import re
from typing import Tuple


class DNA(object):
    """
    Class for holding information related to any input double-stranded DNA sequence used for assembly.

    I wanted to use scikit-bio's DNA class but you can't subclass it... this is a barebone implementation to get
    assembly done. A few things I want done:

    * Circular representation of plasmids
    * Special characters representing overhangs to facilitate golden gate assembly
    * Sequence features that can be inherited by cloning products

    I would use scikit-bio's DNA class for all your analysis and stuff.
    """

    # --- Class variables --- #
    # figure out this compares to static variables in traditional OOP languages
    dna_basepairs = {
        'A': 'T', 
        'T': 'A', 
        'C': 'G', 
        'G': 'C'
    }

    def __init__(self, sequence: str, strand: int = 1, entity_id=None, name=None, features=None, description=None, source=None, metadata=None):
        """
        :param sequence: string representation of dsDNA 5' -> 3'
        :param strand: 1 for sense, -1 for antisense
        :param metadata: dict containing random metadata about the plasmid
        """
        self.entity_id = entity_id
        self.name = name
        self.strand = strand
        self.sequence = sequence
        self.features = features
        self.description = description
        self.source = source
        self.feature_map = None  # Populated by map_features method

    def __repr__(self) -> str:
        """Retrieves a string representation of this DNA object

        Returns:
            str: representation of this DNA object
        """
        return f'DNA:\t{self.entity_id}\t\t{self.name}\t\tLength: {len(self.sequence)}\n' \
               f'Description: {self.description[:50]}{"..." if len(self.sequence) > 50 else ""}\n' \
               f'{self.sequence[:60]}{"..." if len(self.sequence) > 60 else ""}\n'

    def __len__(self) -> int:
        """Retrieves the length of the sequence from this DNA object

        Returns:
            int: length of this DNA sequence
        """
        return len(self.sequence)

    # --- Property Setters --- #

    @property
    def sequence(self) -> str:
        """Retrieves the sequence from this DNA object

        Returns:
            str: sequence of this DNA object
        """
        return self._sequence

    @sequence.setter
    def sequence(self, input_sequence: str) -> None:
        """Assigns a new sequence to this DNA object if the provided sequence
        passes validation.

        Args:
            input_sequence (str): New sequence to assign to this DNA object

        Raises:
            SequenceException: Excepton when a character other than A, T, C, or G
            is present in the sequence string
        """
        
        if re.fullmatch("[ATCG]+", input_sequence.upper()) is not None:
            self._sequence = input_sequence.upper()
        else:
            raise SequenceException("DNA sequences can only contain ATCG.")

    @property
    def strand(self) -> int:
        """Retrieves the strand orientation for the sequence of this DNA object

        Returns:
            int: 1 for sense, -1 for antisense
        """
        return self._strand

    @strand.setter
    def strand(self, strand_int: int) -> None:
        """Assigns a value to indicate strand orientation for the sequence of 
        this DNA object

        Args:
            strand_int (int): 1 for sense, -1 for antisense

        Raises:
            Exception: Error when a value other than 1 or -1 is assigned
        """
        if strand_int in [1, -1]:
            self._strand = strand_int
        else:
            raise Exception('strand attribute can only take the values {-1, 1}')

    # --- Methods --- #
    # should this become property?
    def complement(self) -> str:
        """Generate the complement to the sequence of this DNA object

        Returns:
            str: complement to the sequence of this DNA object
        """

        """
        Pulling code out of list comprehension:

        complement_seq = ""

        for base in self.sequence:
            if base.upper() in self.dna_basepairs.keys():
                complement_seq += self.dna_basepairs.get(base.upper())
            else:
                complement_seq += base
        
        return complement_seq
        """
        return ''.join([self.dna_basepairs.get(base.upper()) if base.upper() in self.dna_basepairs.keys() else base for base in self.sequence])

    # could potentially make this algorithm more effecient
    def reverse_complement(self) -> str:
        """Generate the reverse complement to the sequence of this DNA object

        Returns:
            str: reverse complement to the sequence of this DNA object
        """

        """
        More effecient approach outlined below. Time complexity reduced from O(n^2) to O(n)

        rev_complement_seq = ""

        for base in self.sequence:
            if base.upper() in self.dna_basepairs.keys():
                rev_complement_seq = self.dna_basepairs.get(base.upper()) + rev_complement_seq
            else:
                rev_complement_seq = base + rev_complement_seq
        
        return rev_complement_seq
        """
        return self.complement()[::-1]

    # why is this empty?
    def translate(self):
        """
        Translate DNA sequence into protein
        :return:
        """
        pass

    # incorrectly states in docstring and returns feature map
    def map_features(self) -> None:
        """
        Maps seqeunce features onto the DNA sequence
        :return: {(start, end): feature}
        """
        feature_dict = dict()

        for feature in self.features:
            regex_pattern = f'{feature.sequence}|{feature.reverse_complement()}'
            matches = re.finditer(regex_pattern, self.sequence)

            for match in matches:
                strand = 1 if match.group() == feature.sequence else -1
                feature_dict_key = (match.start(), match.end(), strand)

                # Store features with same locations in list
                if feature_dict_key in feature_dict.keys():
                    feature_dict[feature_dict_key].append(feature)
                else:
                    feature_dict[(match.start(), match.end(), strand)] = [feature]

        self.feature_map = feature_dict

    # make tuple typing more specific
    def find_cut_indicies(self, cut_position, rxn_enzyme) -> Tuple:
        """
        Finds the cut indicies for a rxn_enzyme found at cut_position on sequence
        :param cut_position: cut position index
        :param rxn_enzyme: BioPython Restriction object
        :return:
        """

        # Find restriciton site on coding or non-coding strand
        # todo: make this more permissive to allow fuzzy input cut poisitions
        rxn_match = rxn_enzyme.compsite.search(self.sequence, cut_position, cut_position + len(rxn_enzyme.site))

        if rxn_match is None:
            raise RestrictionSiteException('Restriction site not found!')

        # If rxn_match returns a Match, it has to be rxn_enzyme.site or its reverse complement
        if rxn_match.group(1) == rxn_enzyme.site:
            strand = 1
        else:
            strand = -1

        fst5, fst3, snd5, snd3, rxn_site = rxn_enzyme.charac

        # --- Find overhang if overhang exists --- #

        # Add fst5/3 from regex position for cutsite if rxn site is on coding, else subtract
        fst5 *= strand
        fst3 *= strand

        if strand == -1:
            cut_index_5 = cut_position + len(rxn_site) + fst5
            cut_index_3 = cut_position + fst3
        else:
            cut_index_5 = cut_position + fst5
            cut_index_3 = cut_position + len(rxn_site) + fst3

        return cut_index_5, cut_index_3, strand


class SequenceException(Exception):
    """Custom exception for sequence concent errors"""
    pass


class RestrictionSiteException(Exception):
    """Custom exception for restriction site location errors"""
    pass
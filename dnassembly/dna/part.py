#! /usr/bin/env python3

import re

from .dna import DNA, SequenceException

class Part(DNA):
    """
    Subclass to hold special attributes related to a DNA part. Parts are the individual units of DNA that come together
    to form an assembly and are linear by definition.
    """

    def __init__(self, sequence, id=None, name=None, features=None, description=None, source=None):
        """
        :param name: name of this Part
        :param sequence: string representation of dsDNA 5' -> 3'
        :param overhang_5: tuple of (overhang length, )
        """
        super().__init__(sequence, id=id, name=name, features=features, description=description, source=source)
        self.overhang_5 = None  # (4, True)
        self.overhang_3 = None

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, input_sequence):
        """
        Validate input DNA sequence. Part sequences allow special characters to denote breaks in the top strand
        ▼ (\u25BC) and bottom strand ▲ (\u25B2) of dsDNA to create sticky ends. These symbols are typically seen in
        diagrams outlining restriction sites.

        For example, a typical part in the MoClo system digested with the type II restriction enzyme BsaI:

        5′ - agtcGGTCTCaTATGggttctNNN...NNNgagaggATCCtGAGACCagac - 3'
             ||||||||||||||||||||||||   |||||||||||||||||||||||
        3′ - tcagCCAGAGtATACccaagcNNN...NNNctctccTAGGaCTCTGGtctg - 5'

        Produces the following product when digested:

        5′ -            TATGggttctNNN...NNNgagagg                - 3'
                            |||||||||   ||||||||
        3′ -                ccaagcNNN...NNNctctccTAGG            - 5'

        This collapses into the string representation:

        5′ -           TATG▲ggttctNNN...NNNgagagg▼TAGG           - 3'

        The ▼/▲ symbols make sense if you think of them as denoting which strand of the linear dsDNA gets cut:

        ▼ cuts the (+)-strand / coding strand.
        ▲ cuts the (-)-strand / non-coding strand.

        These symbols will make golden gate part assembly super simple with a little regex and string manipulation.

        :param input_sequence: string representation of input DNA sequence
        :return:
        """

        split_sequence = re.split('[\u25B2\u25BC]', input_sequence)

        # No sticky ends
        if len(split_sequence) == 1:
            self._sequence = input_sequence.upper()

        # Throw exception if more than two nicks are encountered...
        elif len(split_sequence) > 3:
            raise SequenceException('A part can only have two sticky ends.')

        # Make sure sequence is only ATCG between nick markers
        else:
            clean_sequence_parts = [a for a in split_sequence if a is not '']
            if all([re.fullmatch('[ATCG]+', seq.upper()) is not None for seq in clean_sequence_parts]):
                self._sequence = input_sequence.upper()
            else:
                raise SequenceException('DNA sequences can only contain ATCG or nick markers (aka sticky end triangles)')

    @classmethod
    def define_overhangs(cls, sequence, l_overhang_strand=None, l_overhang_bases=0, r_overhang_strand=None, r_overhang_bases=0,
                         id=None, name=None, features=None, description=None, source=None):
        """
        Alternative constructor where sticky ends may be defined in a more conventional/intuitive way... a la BioPython.
        Left and right refer to the 5' and 3' ends of the DNA, respectively. The left and right terminology is used here
        because 5'/3' also refer to the overhang type.

        For the following example part:

        5′ -            TATGggttctNNN...NNNgagagg                - 3'
        Left                |||||||||   ||||||||                Right
        3′ -                ccaagcNNN...NNNctctccTAGG            - 5'

        Constructor inputs:
        * sequence: TATGggttctNNN...NNNgagaggTAGG
        * l_overhang_strand: 5
        * l_overhang_bases: 4
        * r_overhang_strand: 5
        * l_overhang_bases: 4

        This collapses into the string representation:

        5′ -           TATG▲ggttctNNN...NNNgagagg▼TAGG           - 3'

        :param l_overhang_strand: which strand forms the overhang: int({ 3 | 5 }) or None
        :param l_overhang_bases: how many bases form the overhang: int()
        :param r_overhang_strand: which strand forms the overhang: int({ 3 | 5 }) or None
        :param r_overhang_bases: how many bases form the overhang: int()
        :return:
        """

        # --- Checks and warnings --- #
        if l_overhang_bases == 0 and r_overhang_bases == 0:
            print('Warning: both overhangs are of length 0. You should just use the regular Part constructor.')
        if l_overhang_bases <= 0:
            raise PartDefinitionException('l_overhang_bases must be 0 or greater!')
        if r_overhang_bases <= 0:
            raise PartDefinitionException('r_overhang_bases must be 0 or greater!')

        # --- Process sequence into Part syntax --- #

        # Low effort solution... Introduce cuts to left/right side of DNA so that desired overhang type remains
        l_overhang_type = {3: '▲', 5: '▼'}
        r_overhang_type = {3: '▼', 5: '▲'}

        # Right overhang
        right_overhang_sequence = ''
        if r_overhang_strand is None:
            r_overhang_bases = 0  # Ignore user input if r_overhang_strand == None
        else:
            if r_overhang_strand not in [3, 5]:
                raise PartDefinitionException(f'r_overhang_strand can only be int(3), int(5), or None! ({type(r_overhang_strand), r_overhang_strand} was provided)')
            right_overhang_sequence = r_overhang_type[r_overhang_strand] + sequence[-r_overhang_bases:]

        # Left overhang
        left_overhang_sequence = ''
        if l_overhang_strand is None:
            l_overhang_bases = 0  # Ignore user input if l_overhang_strand == None
        else:
            if l_overhang_strand not in [3, 5]:
                raise PartDefinitionException(f'l_overhang_strand can only be int(3), int(5), or None! ({type(l_overhang_strand), l_overhang_strand} was provided)')
            left_overhang_sequence = sequence[:l_overhang_bases] + l_overhang_type[l_overhang_strand]

        # Annotated sequence
        annotated_sequence = left_overhang_sequence + sequence[r_overhang_bases:-l_overhang_bases] + right_overhang_sequence

        return cls(annotated_sequence, id=id, name=name, description=description, features=features, source=source)

    def cut(self, cut_position, rxn_enzyme):
        """
        Cut Part at cut_position with rxn_enzyme to yield two new Parts
        This method manipulates Part.sequence and expects a cut position found on Part.sequence!

        :param cut_position: position on Part sequence to find restriction site for rxn_enzyme
        :param rxn_enzyme: BioPython Restriction object of a Restriction Enzyme
        :return: Part_1, Part_2
        """
        cut_index_5, cut_index_3 = self.find_cut_indicies(cut_position, rxn_enzyme)

        # Handle blunt-end cuts
        if cut_index_5 == cut_index_3:
            new_part_1 = self.sequence[:cut_index_5]
            new_part_2 = self.sequence[cut_index_5:]

        else:
            left_cut = min(cut_index_5, cut_index_3)
            right_cut = max(cut_index_5, cut_index_3)
            sticky_end_sequence = self.sequence[left_cut:right_cut]

            if '▲' in sticky_end_sequence or '▼' in sticky_end_sequence:
                return self, None

            part_1_overhang = '▼' if cut_index_5 > cut_index_3 else '▲'
            part_2_overhang = '▲' if cut_index_5 > cut_index_3 else '▼'

            new_part_1 = self.sequence[:cut_index_5] + part_1_overhang + sticky_end_sequence
            new_part_2 = sticky_end_sequence + part_2_overhang + self.sequence[cut_index_5:]

        return Part(new_part_1, id=self.id, name=f'{self.name}-l_product', features=self.features, source=self.id), \
               Part(new_part_2, id=self.id, name=f'{self.name}-r_product', features=self.features, source=self.id)

    def return_pure_sequence(self):
        """
        :return: part sequence without special symbols, only bases
        """
        return ''.join([base for base in self.sequence if base in 'ATCG'])


class Backbone(Part):
    """
    Subclass to hold special attributes related to a linearized backbone
    """

    def __init__(self):
        pass


class Library(Part):
    """
    Subclass to hold special attributes related to degenerate libraries
    """

    def __init__(self):
        pass


class PartDefinitionException(Exception):
    pass
#! /usr/bin/env python3

import re

from .dna import DNA

class Part(DNA):
    """
    Subclass to hold special attributes related to a DNA part. Parts are the individual units of DNA that come together
    to form an assembly and are linear by definition.
    """

    def __init__(self, sequence, entity_id=None, name=None, features=None, description=None, source=None, overhang_5=None,
                 overhang_3=None):
        """
        :param sequence: string representation of dsDNA 5' -> 3'
        :param entity_id: string, systematic identification of this Part
        :param name: string, user-defined name of this Part
        :param features: list of Feature objects
        :param description: string, describe your part
        :param source: string, entity_id of part source or None
        :param overhang_5: tuple, (overhang length, { 3 | 5 | None }) for 5' end of DNA
        :param overhang_3: tuple, (overhang length, { 3 | 5 | None }) for 3' end of DNA
        """
        super().__init__(sequence, entity_id=entity_id, name=name, features=features, description=description, source=source)
        self.overhang_5 = overhang_5
        self.overhang_3 = overhang_3

    def __repr__(self):
        return f'DNA:\t{self.entity_id}\t\t{self.name}\t\tlength: {len(self.sequence)}\n' \
               f'{self.sequence[:25]} ... {self.sequence[-25:]}\n'

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash((self.sequence, self.overhang_5, self.overhang_3))

    # --- Property Setters --- #

    @property
    def overhang_5(self):
        return self._overhang_5

    @overhang_5.setter
    def overhang_5(self, overhang_5):
        if not (type(overhang_5) is tuple or overhang_5 is None):
            print(overhang_5)
            raise PartDefinitionException('overhang_5 needs to be a tuple or None!')
        if type(overhang_5) is tuple:
            overhang_5_length, overhang_5_strand = overhang_5
            if (overhang_5_length < 0 or overhang_5_length > len(self.sequence)):
                raise PartDefinitionException('overhang_5 length has to be between 0 and the length of your part!')
            if overhang_5_strand not in [3,5,None]:
                raise PartDefinitionException('overhang_5 strand has to be 3, 5, or None!')
        self._overhang_5 = overhang_5

    @property
    def overhang_3(self):
        return self._overhang_3

    @overhang_3.setter
    def overhang_3(self, overhang_3):
        if not (type(overhang_3) is tuple or overhang_3 is None):
            raise PartDefinitionException('overhang_3 needs to be a tuple or None!')
        if type(overhang_3) is tuple:
            overhang_3_length, overhang_3_strand = overhang_3
            if (overhang_3_length < 0 or overhang_3_length > len(self.sequence)):
                raise PartDefinitionException('overhang_3 length has to be between 0 and the length of your part!')
            if overhang_3_strand not in [3,5,None]:
                raise PartDefinitionException('overhang_3 strand has to be 3, 5, or None!')
        self._overhang_3 = overhang_3

    @classmethod
    def define_overhangs(cls, sequence, l_overhang_strand=None, l_overhang_bases=0, r_overhang_strand=None, r_overhang_bases=0,
                         entity_id=None, name=None, features=None, description=None, source=None):
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

        # --- Process overhangs if they exist --- #

        # 5' overhang
        if l_overhang_bases == 0:
            overhang_5 = None
        else:
            if l_overhang_strand not in [3, 5]:
                raise PartDefinitionException(f'l_overhang_strand can only be int(3) or int(5)! ({type(l_overhang_strand), l_overhang_strand} was provided)')
            overhang_5 = (l_overhang_bases, l_overhang_strand)

        # 3' overhand
        if l_overhang_bases == 0:
            overhang_3 = None
        else:
            if r_overhang_strand not in [3, 5]:
                raise PartDefinitionException(f'r_overhang_strand can only be int(3) or int(5)! ({type(r_overhang_strand), r_overhang_strand} was provided)')
            overhang_3 = (r_overhang_bases, r_overhang_strand)

        return cls(sequence, entity_id=entity_id, name=name, description=description, features=features, source=source, overhang_3=overhang_3, overhang_5=overhang_5)

    # --- Methods --- #

    def cut(self, cut_position, rxn_enzyme):
        """
        Cut Part at cut_position with rxn_enzyme to yield two new Parts. This method manipulates Part.sequence and
        expects a cut position found on Part.sequence!

        :param cut_position: position on Part sequence to find restriction site for rxn_enzyme
        :param rxn_enzyme: BioPython Restriction object of a Restriction Enzyme
        :return: Part_1, Part_2
        """
        cut_index_5, cut_index_3, strand = self.find_cut_indicies(cut_position, rxn_enzyme)

        # Check if restriction site has already been processed
        stick_end_offset = 0 if self.overhang_3 is None else self.overhang_3[0]
        if max(cut_index_5, cut_index_3) >= (len(self.sequence) - stick_end_offset):
            print('Restriction site already processed!')
            return self, None

        # Handle blunt-end cuts
        if cut_index_5 == cut_index_3:
            new_part_1 = self.sequence[:cut_index_5]
            new_part_2 = self.sequence[cut_index_5:]
            part_1_overhang_3 = None
            part_2_overhang_5 = None

        # Handle sticky ends
        else:
            overhang = 5 if (cut_index_3 - cut_index_5) * strand > 0 else 3
            part_1_overhang_3 = part_2_overhang_5 = (abs(cut_index_5 - cut_index_3), overhang)

            new_part_1 = self.sequence[:max(cut_index_3, cut_index_5)]
            new_part_2 = self.sequence[min(cut_index_3, cut_index_5):]

        return Part(new_part_1, entity_id=self.entity_id, name=f'{self.name}-l_product', features=self.features, source=self.entity_id,
                    overhang_5=self.overhang_5, overhang_3=part_1_overhang_3),\
               Part(new_part_2, entity_id=self.entity_id, name=f'{self.name}-r_product', features=self.features, source=self.entity_id,
                    overhang_5=part_2_overhang_5, overhang_3=self.overhang_3)


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
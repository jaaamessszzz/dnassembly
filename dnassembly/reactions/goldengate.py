#! /usr/bin/env python3

import re

from .cloning import StickyEndAssembly, ReactionDefinitionException
from ..dna import Plasmid
from ..utils import pairwise

class GoldenGate(StickyEndAssembly):
    """
    Perfoms golden gate assembly using input parts

    :return:
    """

    def __init__(self, input_dna_list, restriction_enzyme_list):
        super(GoldenGate, self).__init__(input_dna_list, restriction_enzyme_list)
        self.verify_parts()

    # --- Setters --- #

    @property
    def restriction_enzyme_list(self):
        return self._restriction_enzyme_list

    @restriction_enzyme_list.setter
    def restriction_enzyme_list(self, enzyme):
        if len(enzyme) > 1:
            raise ReactionDefinitionException('There can only be one enzyme in a Golden Gate reaction!')
        elif (len(enzyme.site) + enzyme.fst3 < len(enzyme.site) or enzyme.fst5 <len(enzyme.site)):
            raise ReactionDefinitionException('Type II Restriction Enzymes are required for Golden Gate Reactions!')
        else:
            self._restriction_enzyme_list = [type.value]

    # --- Methods --- #

    def verify_parts(self):
        """
        Verify input plasmids only contain two restriction sites of the proper type and orientation
        Find restriction sites in sequence, then deque full 180 and find restriction sits again
        :return:
        """
        # todo: finsih writing this...
        for dna_entity in self.input_dna_list:
            if isinstance(dna_entity, Plasmid):
                for rxn_enzyme in self.restriction_enzyme_list:

                    cutsites = list()

                    # Finds all non-overlapping restriction sites... so don't overlap restriction sites
                    for cutsite in re.finditer(rxn_enzyme.compsite, dna_entity.sequence):
                        cut_position = cutsite.start(0)

                        if dna_entity.sequence[cut_position:cut_position + len(rxn_enzyme.site)] == rxn_enzyme.site:
                            strand = 1
                        else:
                            strand = -1
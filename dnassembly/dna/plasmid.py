#! /usr/bin/env python3
import re

from .dna import DNA
from .part import Part
from ..utils import cycle_in_frames

class Plasmid(DNA):
    """
    Completed product of a DNA assembly. Plasmids are circular by definition, an all methods of Plasmid take this into
    account (e.g. looking for features)

    Look into
    from collections import deque
    for circular plasmid representation
    """

    def __init__(self, sequence, id=None, name=None, features=None, description=None, source=None):
        super().__init__(sequence, id=id, name=name, features=features, description=description, source=source)

    def __repr__(self):
        return f'\nPlasmid:\t{self.id}\t\t{self.name}\t\tlength: {len(self.sequence)}\n' \
               f'{self.description:<80}\n' \
               f'{self.sequence[:25]}...\n'

    def linearize(self, cut_position, rxn_enzyme):
        """
        Linearize Plasmid at a given cut site with rxn_enzyme to produce a Part
        :return:
        """

        cut_index_5, cut_index_3 = self.find_cut_indicies(cut_position, rxn_enzyme)

        # --- Deque so cut first blunt/sticky end is at end of linear sequence --- #

        sequence_intermediate = self.sequence
        offset = max(cut_index_5, cut_index_3) - len(sequence_intermediate)
        sequence_intermediate = sequence_intermediate[offset:] + sequence_intermediate[:offset]

        # If sticky end exists, add repeat sticky end to beginning and figure out which strand forms the overhang
        sticky_end_size = abs(cut_index_5 - cut_index_3)
        sequence_intermediate = sequence_intermediate[-sticky_end_size:] + sequence_intermediate

        if cut_index_5 > cut_index_3:
            overhang_strand = 5
        elif cut_index_5 < cut_index_3:
            overhang_strand = 3
        else:
            overhang_strand = None

        input_dna = Part.define_overhangs(sequence_intermediate,
                                          l_overhang_strand=overhang_strand,
                                          l_overhang_bases=sticky_end_size,
                                          r_overhang_strand=overhang_strand,
                                          r_overhang_bases=sticky_end_size,
                                          id=self.id,
                                          name=f'{self.name}_linear',
                                          features=self.features,
                                          source=self.id)

        return input_dna

    def digest_plasmid(self, rxn_enzyme_list):
        """
        Digests plasmid into parts
        :return:
        """


    def map_features(self):
        """
        Maps seqeunce features onto the plasmid sequence. This overrides the map_feature method in DNA to ensure
        features that span the beginning/end of the sequence are annotated.

        :return: {(start, end): feature}
        """
        pass


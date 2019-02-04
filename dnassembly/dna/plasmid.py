#! /usr/bin/env python3

import re

from .dna import DNA
from .part import Part


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
        self.feature_map = None  # Populated by map_features method

    def __repr__(self):
        return f'\nPlasmid:\t{self.id}\t\t{self.name}\t\tlength: {len(self.sequence)}\n' \
               f'{self.description:<80}\n' \
               f'{self.sequence[:25]}...\n'

    # --- Methods --- #

    def linearize(self, cut_position, rxn_enzyme):
        """
        Linearize Plasmid at a given cut site with rxn_enzyme to produce a Part
        :return:
        """

        cut_index_5, cut_index_3 = self.find_cut_indicies(cut_position, rxn_enzyme)

        # --- Deque so cut first blunt/sticky end is at end of linear sequence --- #

        sequence_intermediate = self.sequence
        offset = max(cut_index_5, cut_index_3)
        sequence_intermediate = sequence_intermediate[offset:] + sequence_intermediate[:offset]

        # If sticky end exists, add repeat sticky end to beginning and figure out which strand forms the overhang
        sticky_end_size = abs(cut_index_5 - cut_index_3)
        sequence_intermediate = sequence_intermediate[-sticky_end_size:] + sequence_intermediate

        if cut_index_5 == cut_index_3:
            overhang_strand = None
        else:
            overhang_strand = 5 if cut_index_5 < cut_index_3 else 3

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

    def map_features(self):
        """
        Maps seqeunce features onto the plasmid sequence. This overrides the map_feature method in DNA to ensure
        features that span the beginning/end of the sequence are annotated.

        :return: {(start, end, strand): feature} where strand = 1 for 5' -> 3', strand = -1 for 3' -> 5'
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





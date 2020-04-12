#! /usr/bin/env python3

import re

from .dna import SequenceException
from .part import Part


class Feature(Part):
    """
    A unique sequence of DNA that performs some biological function of interest
    """
    def __init__(self, sequence, feature_type, strand, entity_id=None, name=None, features=None, description=None,
                 forward_color=None, reverse_color=None):
        super(Feature, self).__init__(sequence, entity_id=entity_id, name=name, features=features, description=description)
        self.feature_type = feature_type
        self.strand = strand
        self.forward_color = forward_color
        self.reverse_color = reverse_color

    def __repr__(self):
        return f'Feature:\t{self.entity_id}\t\t{self.name}\t\tlength: ({len(self.sequence)})\n{self.sequence[:25]}...\n'

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash((self.sequence, self.feature_type, self.entity_id))

    # --- Property Setters --- #

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, input_sequence):
        """
        Validate input DNA sequence. Features can have gaps or degenerate recognition sites.

        :param input_sequence: string representation of input DNA sequence
        :return:
        """
        # Stupidity check
        if re.fullmatch('[.]+', input_sequence.upper()) is not None:
            raise SequenceException('Why would you do that...')

        if re.fullmatch('[ATCG\.RYSWKMBDHVN]+', input_sequence.upper()) is not None:
            self._sequence = input_sequence.upper()
        else:
            raise SequenceException('Features can only contain valid IUPAC nucleotide symbols or "." to denote gaps')
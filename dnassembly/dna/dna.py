#! /usr/bin/env python3

import re

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

    dna_basepairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence, entity_id=None, name=None, features=None, description=None, source=None):
        """
        :param sequence: string representation of dsDNA 5' -> 3'
        """
        self.entity_id = entity_id
        self.name = name
        self.sequence = sequence
        self.features = features
        self.description = description
        self.source = source
        self.feature_map = None  # Populated by map_features method

    def __repr__(self):
        return f'DNA:\t{self.entity_id}\t\t{self.name}\t\tlength: {len(self.sequence)}\n' \
               f'{self.sequence[:25]} ... {self.sequence[-25:]}\n'

    # --- Property Setters --- #

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, input_sequence):
        """
        Validate input DNA sequence

        :param input_sequence: string representation of input DNA sequence
        :return:
        """
        if re.fullmatch('[ATCG]+', input_sequence.upper()) is not None:
            self._sequence = input_sequence.upper()
        else:
            raise SequenceException('DNA sequences can only contain ATCG.')

    # --- Methods --- #

    def complement(self):
        """
        Generate the complement to the DNA sequence
        :return:
        """
        return ''.join([self.dna_basepairs.get(base.upper()) if base.upper() in self.dna_basepairs.keys() else base for base in self.sequence])

    def reverse_complement(self):
        """
        Generate the reverse complement to the DNA sequence
        :return:
        """
        return self.complement()[::-1]

    def translate(self):
        """
        Translate DNA sequence into protein
        :return:
        """
        pass

    def map_features(self):
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

    def find_cut_indicies(self, cut_position, rxn_enzyme):
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


class Feature(DNA):
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


class SequenceException(Exception):
    pass


class RestrictionSiteException(Exception):
    pass
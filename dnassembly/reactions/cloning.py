#! /usr/bin/env python3

import re
from pprint import pprint
from collections import OrderedDict

from Bio.Restriction.Restriction import RestrictionType

from ..dna import DNA, Plasmid, Part
from ..utils import pairwise


class CloningReaction(object):
    """
    An instance of CloningReaction represents a single-step DNA manipulation that you would carry out in an Eppendorf or
    PCR strip tube (you should really be doing everything in strip tubes though). Examples of a cloning reaction would
    be a test restriction digest or Golden Gate reaction.

    Assembly reactions such as Gibson and GoldenGate subclass CloningReaction to impose their own biological constraints.
    """

    def __init__(self, input_dna_list, restriction_enzyme_list):
        """
        :param input_dna_list: list of DNA entities for a particular cloning reaction
        """
        self.input_dna_list = input_dna_list
        self.restriction_enzyme_list = restriction_enzyme_list
        self.digest_pool = None  # Populated by self.digest()

    @property
    def input_dna_list(self):
        return self._input_dna_list

    @input_dna_list.setter
    def input_dna_list(self, input_dna_list):
        """
        Adds DNA entities to this reaction
        :param dna:
        :return:
        """
        if type(input_dna_list) is not list:
            input_dna_list = [input_dna_list]

        dna_list = list()
        for dna in input_dna_list:
            if isinstance(dna, DNA):
                dna_list.append(dna)
            else:
                raise DNAException('Only DNA can be added as substrates to cloning reactions!')
        self._input_dna_list = dna_list

    @property
    def restriction_enzyme_list(self):
        return self._restriction_enzyme_list

    @restriction_enzyme_list.setter
    def restriction_enzyme_list(self, input_enzyme_list):
        if type(input_enzyme_list) is not list:
            input_enzyme_list = [input_enzyme_list]

        enzyme_list = list()
        for enzyme in input_enzyme_list:
            if isinstance(enzyme, RestrictionType):
                enzyme_list.append(enzyme)
            else:
                raise EnzymeException('Only BioPython Restriction Enzymes can be added to catalyze reactions!')
        self._restriction_enzyme_list = enzyme_list

    def find_restriction_sites(self, input_dna):
        """
        Find restriction sites in input_dna for enzymes in restriction_enzyme_list
        :param input_dna: DNA entity
        :return:
        """
        # Store enzyme cut information so we can process all cuts from all enzymes on full sequence in serial
        rxn_enzyme_cuts = dict()

        # Identify restriction sites for all enzymes in the reaction
        for rxn_enzyme in self.restriction_enzyme_list:

            # Finds all non-overlapping restriction sites... so don't overlap restriction sites
            for cutsite in re.finditer(rxn_enzyme.compsite, input_dna.sequence):

                cut_position = cutsite.start(0)
                if input_dna.sequence[cut_position:cut_position + len(rxn_enzyme.site)] == rxn_enzyme.site:
                    strand = 1
                else:
                    strand = -1

                # Store enzyme cut information
                rxn_enzyme_cuts[cutsite.start(0)] = {'enzyme': rxn_enzyme,
                                                     'strand': strand}

        return rxn_enzyme_cuts

    def digest_recursively(self):
        """
        Digest all DNA entities in input_dna_list with with all enzymes in restriction_enzyme_list to produce Parts.

        This method will find all restriction sites on a part and recursively break the part down (3' -> 5') until all
        restriction sites are processed to produce all digestion products. Plasmids are linearized into Parts before
        digestion.

        Why do this recursively? Most edge cases are covered without additional bs checks.
        * Already processed restriction sites caught early-on
        * Linearizing plasmids before digestion removes issues with restriction sites spanning start/end of sequence

        :return: list of Parts() keyed to their source DNA
        """
        digestion_products = list()

        # For each input DNA entity
        for input_dna in self.input_dna_list:

            # --- Find restriction sites in input_dna --- #

            rxn_enzyme_cuts = self.find_restriction_sites(input_dna)

            # --- Break input_dna into Parts --- #

            # Linearize Plasmid if necessary
            if isinstance(input_dna, Plasmid):

                # Pick a random site to linearize Plasmid
                cut_position = next(iter(rxn_enzyme_cuts.keys()))
                input_dna = input_dna.linearize(cut_position, rxn_enzyme_cuts[cut_position]['enzyme'])

                # Redo restriction site search to find any restriction sites that may have straddled start/end
                # Efficient? No. Easy and avoids problems? Yes.
                rxn_enzyme_cuts = self.find_restriction_sites(input_dna)

            # Create OrderedDict sorted by cut site order (Recipe from OrderedDict documentation!)
            sorted_rxn_enzyme_cuts = OrderedDict(reversed(sorted(rxn_enzyme_cuts.items(), key=lambda t: t[0])))

            # Proceed through sequence in reverse and break down into Parts
            source_plasmid = input_dna.id
            for cutsite_position, cutsite_info in sorted_rxn_enzyme_cuts.items():
                input_dna, part_temp = input_dna.cut(cutsite_position, cutsite_info['enzyme'])

                # Part.cut returns part_temp as None if the current restriction site has already been processed
                if part_temp is None:
                    continue

                part_temp.source = source_plasmid
                digestion_products.append(part_temp)

            # input_dna is the final part after all cuts have been processed
            input_dna.source = source_plasmid
            digestion_products.append(input_dna)

        # Populate digest pool
        self.digest_pool = digestion_products

    def digest_by_slicing(self):
        """
        Iterate through rxn_enzyme_cuts in order and pull sequences from DNA to produce new Parts.
        Super dirty initial implementation, we'll clean that up later (he said 5 years ago)...
        There are a lot of edge cases I have to check using this method, specifically catching restriction sites that
        span the start/end of a plasmid sequence...
        and if a plasmid only has one cut...
        :return:
        """
        def process_5(cut_index_5, cut_index_3):
            cut_index_difference = cut_index_5 - cut_index_3
            if cut_index_difference == 0:
                return None
            else:
                overhang_strand = 3 if cut_index_difference > 0 else 5
                return (abs(cut_index_difference), overhang_strand)

        def process_3(cut_index_5, cut_index_3):
            cut_index_difference = cut_index_5 - cut_index_3
            if cut_index_difference == 0:
                return None
            else:
                overhang_strand = 3 if cut_index_difference < 0 else 5
                return (abs(cut_index_difference), overhang_strand)

        # Nested because this function really should not be exposed, defeats the point of CloningReaction objects
        def make_part(template_dna, left_cut, right_cut, plasmid_span=False):
            """
            Use {cut_position: rxn_enzyme_dict} to create a new part from input_dna
            :param left_cut: {cut_position: {'enzyme': BioPython Restriction, 'strand': {1 | -1} } }
            :param right_cut: {cut_position: {'enzyme': BioPython Restriction, 'strand': {1 | -1} } }
            :param plasmid_span: assemble
            """

            "There are subtle differences in how 5' and 3' ends are processed which makes generalizing code difficult."
            # --- Process 5' of new Part --- #

            if left_cut[0] == 'Start':
                part_overhang_5 = left_cut[1]
                left_cut_index = 0
                left_cut_indicies = (0,0)

            else:
                cut_position_left, rxn_enzyme_dict_left = left_cut
                cut_index_5, cut_index_3 = template_dna.find_cut_indicies(cut_position_left, rxn_enzyme_dict_left['enzyme'])
                left_cut_indicies = (cut_index_5, cut_index_3)
                left_cut_index = min(left_cut_indicies)

                if isinstance(template_dna, Plasmid):
                    # Process restriction site
                    part_overhang_5 = process_5(cut_index_5, cut_index_3)
                else:
                    # Handle already-processed restriction sites at part terminals (e.g. digested BsaI part backbone)
                    if min(cut_index_5, cut_index_3) <= (template_dna.overhang_5[0]):
                        part_overhang_5 = template_dna.overhang_5
                    # Process restriction site
                    else:
                        part_overhang_5 = process_5(cut_index_5, cut_index_3)

            # --- Process 3' of new part --- #

            if right_cut[0] == 'End':
                part_overhang_3 = right_cut[1]
                right_cut_index = len(template_dna.sequence)
                right_cut_indices = (right_cut_index, right_cut_index)
            else:
                cut_position_right, rxn_enzyme_dict_right = right_cut
                cut_index_5, cut_index_3 = template_dna.find_cut_indicies(cut_position_right, rxn_enzyme_dict_right['enzyme'])
                right_cut_indices = (cut_index_5, cut_index_3)
                right_cut_index = max(right_cut_indices)

                if isinstance(template_dna, Plasmid):
                    # Process restriction site
                    part_overhang_3 = process_3(cut_index_5, cut_index_3)
                else:
                    # Handle already-processed restriction sites at part terminals (e.g. digested BsaI part backbone)
                    if max(cut_index_5, cut_index_3) >= (len(template_dna.sequence) - template_dna.overhang_3[0]):
                        part_overhang_3 = template_dna.overhang_3
                    # Process restriction site
                    else:
                        part_overhang_3 = process_3(cut_index_5, cut_index_3)

            # --- Check for dud parts --- #
            "These pop up as a result of checking for parts across plasmid sequence end/start."

            # Match between restriction sites directly adjacent to previously processed start/end of part
            if not isinstance(template_dna, Plasmid) and max(left_cut_indicies) >= min(right_cut_indices):
                return None
            # Duplicate part
            if min(left_cut_indicies) == 0 and max(right_cut_indices) == len(template_dna.sequence):
                return None

            # --- Make new Part --- #
            if plasmid_span:
                new_part_sequence = template_dna.sequence[left_cut_index:] + template_dna.sequence[:right_cut_index]
            else:
                new_part_sequence = template_dna.sequence[left_cut_index:right_cut_index]
            new_part = Part(new_part_sequence, id=input_dna.id, name=input_dna.name, features=input_dna.features,
                            description=input_dna.name, source=input_dna.id, overhang_5=part_overhang_5,
                            overhang_3=part_overhang_3)
            print('I made dis.', left_cut_index, right_cut_index)
            return new_part

        digestion_products = list()

        # For each input DNA entity
        for input_dna in self.input_dna_list:

            # --- Find restriction sites in input_dna --- #

            rxn_enzyme_cuts = self.find_restriction_sites(input_dna)
            sorted_rxn_enzyme_cuts = OrderedDict(sorted(rxn_enzyme_cuts.items(), key=lambda t: t[0]))

            # If Part/DNA, add overhang information to start/end of OrderedDict
            if not isinstance(input_dna, Plasmid):
                # Start (5')
                sorted_rxn_enzyme_cuts.update({'Start': input_dna.overhang_5})
                sorted_rxn_enzyme_cuts.move_to_end('Start', last=False)
                # End (3')
                sorted_rxn_enzyme_cuts.update({'End': input_dna.overhang_3})
                sorted_rxn_enzyme_cuts.move_to_end('End', last=True)

            # Duplicate cut if there's only one site in a Plasmid
            if isinstance(input_dna, Plasmid) and len(sorted_rxn_enzyme_cuts) == 1:
                sorted_rxn_enzyme_cuts.update(rxn_enzyme_cuts)

            # --- Perform assemblies --- #

            for left_cut, right_cut in pairwise(sorted_rxn_enzyme_cuts.items()):
                new_part = make_part(input_dna, left_cut, right_cut)
                digestion_products.append(new_part)

            # Perform last part assembly across sequence start/end for Plasmids
            if isinstance(input_dna, Plasmid):

                dict_keys = [a for a in sorted_rxn_enzyme_cuts.keys()]
                right_cut = (dict_keys[0], sorted_rxn_enzyme_cuts[dict_keys[0]])
                left_cut = (dict_keys[-1], sorted_rxn_enzyme_cuts[dict_keys[-1]])
                plasmid_span_part = make_part(input_dna, left_cut, right_cut, plasmid_span=True)

                # Check restriction sites against last part to check for sites that may have spanned gap
                rxn_enzyme_cuts = self.find_restriction_sites(plasmid_span_part)

                if not rxn_enzyme_cuts:
                    digestion_products.append(plasmid_span_part)

                else:
                    sorted_rxn_enzyme_cuts = OrderedDict(sorted(rxn_enzyme_cuts.items(), key=lambda t: t[0]))

                    # Start (5')
                    sorted_rxn_enzyme_cuts.update({'Start': plasmid_span_part.overhang_5})
                    sorted_rxn_enzyme_cuts.move_to_end('Start', last=False)
                    # End (3')
                    sorted_rxn_enzyme_cuts.update({'End': plasmid_span_part.overhang_3})
                    sorted_rxn_enzyme_cuts.move_to_end('End', last=True)

                    for left_cut, right_cut in pairwise(sorted_rxn_enzyme_cuts.items()):

                        plasmid_span_parts = make_part(plasmid_span_part, left_cut, right_cut)
                        digestion_products.append(plasmid_span_parts)

        self.digest_pool = [a for a in digestion_products if a is not None]

    def interpret_biopython_rxn_enzymes(self, rxn_enzyme):
        """
        Pro: BioPython provides classes for individual enzymes from REBASE
        Con: BioPython's Restriction class has not been applied outside of analysis

        We are here to fix that.

        From what I have been able to gather, non-obvious properties of a given enzyme:
        * compsite: regex for restriction site
        * fst3:     position of first cut on 3' strand counting back from end of restriction site
        * fst5:     position of first cut on 5' strand counting forward from beginning of restriciton site
        * charac:   (fst5, fst3, snd5, snd3, compsite)

        :return:
        """
        pass


class DNAException(Exception):
    pass


class EnzymeException(Exception):
    pass
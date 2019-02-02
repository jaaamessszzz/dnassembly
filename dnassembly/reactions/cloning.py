#! /usr/bin/env python3

import re
from pprint import pprint
from collections import OrderedDict

from Bio.Restriction.Restriction import RestrictionType

from ..dna import DNA, Plasmid
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

    def digest(self):
        """
        Digest all DNA entities in input_dna_list with with all enzymes in restriction_enzyme_list to produce Parts.


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
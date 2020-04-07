#! /usr/bin/env python3

import re
from itertools import combinations, product
from collections import OrderedDict

import networkx
from Bio.Restriction.Restriction import RestrictionType

from ..dna import DNA, Plasmid, Part
from dnassembly.utils.utils import pairwise


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
        :param restriction_enzyme_list: list of BioPython Restriction Enzymes for a particular cloning reaction
        :param _digested: has a digestion been run on input_dna_list?
        """
        self.input_dna_list = input_dna_list
        self.restriction_enzyme_list = restriction_enzyme_list
        self._digested = False  # Set to True is successfully _digested by self.digest()

    # --- Property Setters --- #

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
                raise ReactionDefinitionException('Only DNA can be added as substrates to cloning reactions!')
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
                raise ReactionDefinitionException('Only BioPython Restriction Enzymes can be added to catalyze reactions!')
        self._restriction_enzyme_list = enzyme_list

    # --- Methods --- #

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

    def digest(self, recursive=True):
        """
        Digest DNA in input_dna_list using slicing (default) or the recursive method.
        :param recursive: use recursive method for digestion if True
        :return:
        """
        if self._digested is False:
            self._digest_recursively() if recursive else self._digest_by_slicing()
        else:
            raise ReactionDefinitionException('The DNA in this reaction has already been digested!')

    def _digest_recursively(self):
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

            if len(rxn_enzyme_cuts) == 0 and isinstance(input_dna, Plasmid):
                raise AssemblyException(f'Plasmid {input_dna.entity_id} cannot be cut by the restriction enzymes in this reaction!\n'
                                        f'Enzymes: {" ".join([a.__name__ for a in self.restriction_enzyme_list])}')

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
            source_plasmid = input_dna.entity_id
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

        # Replace input_dna_list with new digestion products
        self.input_dna_list = digestion_products
        self._digested = True

    def _digest_by_slicing(self):
        """
        Iterate through rxn_enzyme_cuts in order and pull sequences from DNA to produce new Parts.
        Super dirty initial implementation, we'll clean that up later (he said 5 years ago)...
        There are a lot of edge cases I have to check using this method, specifically catching restriction sites that
        span the start/end of a plasmid sequence...
        and if a plasmid only has one cut...
        :return:
        """
        def process_sticky_end(cut_index_5, cut_index_3, strand):
            cut_index_difference = (cut_index_3 - cut_index_5) * strand
            if cut_index_difference == 0:
                return None
            else:
                overhang_strand = 5 if cut_index_difference > 0 else 3
                return abs(cut_index_difference), overhang_strand

        # Nested because this function really should not be exposed, defeats the point of CloningReaction objects
        def make_part(template_dna, left_cut, right_cut, plasmid_span=False):
            """
            Use {cut_position: rxn_enzyme_dict} to create a new part from input_dna
            There are subtle differences in how 5' and 3' ends are processed which makes generalizing code difficult...

            :param left_cut: {cut_position: {'enzyme': BioPython Restriction, 'strand': {1 | -1} } }
            :param right_cut: {cut_position: {'enzyme': BioPython Restriction, 'strand': {1 | -1} } }
            :param plasmid_span: assemble
            """

            # --- Process 5' of new Part --- #

            if left_cut[0] == 'Start':
                part_overhang_5 = left_cut[1]
                left_cut_index = 0
                left_cut_indicies = (0,0)

            else:
                cut_position_left, rxn_enzyme_dict_left = left_cut
                cut_index_5, cut_index_3, strand = template_dna.find_cut_indicies(cut_position_left, rxn_enzyme_dict_left['enzyme'])
                left_cut_indicies = (cut_index_5, cut_index_3)
                left_cut_index = min(left_cut_indicies)

                if isinstance(template_dna, Plasmid):
                    # Process restriction site
                    part_overhang_5 = process_sticky_end(cut_index_5, cut_index_3, strand)
                else:
                    # Handle already-processed restriction sites at part terminals (e.g. digested BsaI part backbone)
                    stick_end_offset = 0 if template_dna.overhang_5 is None else template_dna.overhang_5[0]
                    if min(cut_index_5, cut_index_3) <= stick_end_offset:
                        part_overhang_5 = template_dna.overhang_5
                    # Process restriction site
                    else:
                        part_overhang_5 = process_sticky_end(cut_index_5, cut_index_3, strand)

            # --- Process 3' of new part --- #

            if right_cut[0] == 'End':
                part_overhang_3 = right_cut[1]
                right_cut_index = len(template_dna.sequence)
                right_cut_indices = (right_cut_index, right_cut_index)
            else:
                cut_position_right, rxn_enzyme_dict_right = right_cut
                cut_index_5, cut_index_3, strand = template_dna.find_cut_indicies(cut_position_right, rxn_enzyme_dict_right['enzyme'])
                right_cut_indices = (cut_index_5, cut_index_3)
                right_cut_index = max(right_cut_indices)

                if isinstance(template_dna, Plasmid):
                    # Process restriction site
                    part_overhang_3 = process_sticky_end(cut_index_5, cut_index_3, strand)
                else:
                    # Handle already-processed restriction sites at part terminals (e.g. digested BsaI part backbone)
                    stick_end_offset = len(template_dna.sequence) if template_dna.overhang_5 is None else (len(template_dna.sequence) - template_dna.overhang_3[0])
                    if max(cut_index_5, cut_index_3) >= stick_end_offset:
                        part_overhang_3 = template_dna.overhang_3
                    # Process restriction site
                    else:
                        part_overhang_3 = process_sticky_end(cut_index_5, cut_index_3, strand)

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
            new_part = Part(new_part_sequence, entity_id=input_dna.entity_id, name=input_dna.name, features=input_dna.features,
                            description=input_dna.name, source=input_dna.entity_id, overhang_5=part_overhang_5,
                            overhang_3=part_overhang_3)
            return new_part

        digestion_products = list()

        # For each input DNA entity
        for input_dna in self.input_dna_list:

            # --- Find restriction sites in input_dna --- #

            rxn_enzyme_cuts = self.find_restriction_sites(input_dna)
            if len(rxn_enzyme_cuts) == 0 and isinstance(input_dna, Plasmid):
                raise AssemblyException(f'Plasmid {input_dna.entity_id} cannot be cut by the restriction enzymes in this reaction!\n'
                                        f'Enzymes: {" ".join([a.__name__ for a in self.restriction_enzyme_list])}')

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

                # Check for restriction sites that may have spanned gap
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

        self.input_dna_list = [a for a in digestion_products if a is not None]
        self._digested = True


class StickyEndAssembly(CloningReaction):
    """
    Perform assemblies using sticky-end methods.

    * Golden Gate
    * BioBrick
    """

    def __init__(self, input_dna_list, restriction_enzyme_list):
        super().__init__(input_dna_list, restriction_enzyme_list)

    def perform_assembly(self, plasmids_only=True, new_id='New_assembly', new_description='New assembly'):
        """
        Performs an assembly provided a pool of Parts.

        Inspiration from pydna to use networkx to represent assemblies as a graph. However I think we can be a lot more
        efficient and enforce our biological constraints at the same time.
        :param plasmids_only: only report cyclic assemblies that fulfill reaction constraints
        :return:
        """

        # --- Sanity Checks --- #

        """
        All DNA entities in input_dna_list should be parts. Any undigested Plasmids with the same resistance marker as
        your assembly product will totally mess up your day.
        """
        if any([isinstance(dna, Plasmid) for dna in self.input_dna_list]):
            raise ReactionDefinitionException('Undigested plasmids cannot be part of an assembly!')

        # --- Assemble Graph --- #

        directed_graph = networkx.MultiDiGraph()

        # Add nodes and edges for each part in digest_pool
        for part in self.input_dna_list:

            sticky_match_l = part.overhang_5
            sticky_match_r = part.overhang_3

            # Process sticky ends into nodes
            # (overhang_sequence, overhang_strand)
            l_node = 'Blunt_Left' if sticky_match_l is None else (part.sequence[:sticky_match_l[0]], sticky_match_l[1])
            r_node = 'Blunt_Right' if sticky_match_r is None else (part.sequence[-sticky_match_r[0]:], sticky_match_r[1])

            directed_graph.add_node(l_node)
            directed_graph.add_node(r_node)

            # Add directed edge between nodes (5' -> 3') and add Part as edge attribute
            directed_graph.add_edge(l_node, r_node, part=part)

        # --- Traverse Graph --- #

        # First: simple_cycles
        graph_cycles = networkx.algorithms.cycles.simple_cycles(directed_graph)
        processed_cycles = list()
        all_possible_assemblies = list()

        for cycle in graph_cycles:

            if not any([set(cycle) == processed_cycle for processed_cycle in processed_cycles]):
                sequences_list = list()

                # Iterate through edges and get sequence up to 3' sticky end
                for node_l, node_r in pairwise(cycle):
                    current_edge = directed_graph[node_l][node_r]
                    new_sequence_list = [current_edge[index]['part'] for index, edge in enumerate(current_edge)]
                    sequences_list.append(new_sequence_list)

                # Get last part of cycle
                last_edge = directed_graph[cycle[-1]][cycle[0]]
                new_sequence_list = [last_edge[index]['part'] for index, edge in enumerate(last_edge)]
                sequences_list.append(new_sequence_list)
                all_possible_assemblies += product(*sequences_list)

        # Actually do assemblies
        intermediate_assemblies = list()
        for assembly in all_possible_assemblies:
            assembly_dict = dict()
            assembly_dict['sequence'] = ''
            assembly_dict['features'] = list()
            assembly_dict['sources'] = list()
            assembly_dict['description'] = list()
            for part in assembly:
                right_sticky_end = part.overhang_3[0]
                assembly_dict['sequence'] += part.sequence[:-right_sticky_end]
                if part.features:
                    assembly_dict['features'] += part.features
                if part.description:
                    assembly_dict['description'].append(part.description)
                assembly_dict['sources'].append(part.source)

            if assembly_dict['sequence'] != '':
                intermediate_assemblies.append(assembly_dict)

        # In plasmids_only=False: all_simple_paths, iterate over all combinations of nodes
        # todo: write this

        # --- Final digestion --- #

        # Omit assemblies that still possess restriction sites for enzymes in restriction_enzyme_list
        complete_assemblies = list()

        for assembly in intermediate_assemblies:
            rxn_site_found = False

            for enzyme in self.restriction_enzyme_list:
                if len(enzyme.compsite.findall(assembly['sequence'])) > 0:
                    rxn_site_found = True

            if rxn_site_found is True:
                continue
            else:
                complete_assemblies.append(assembly)

        # --- Raise exceptions if something went wrong --- #

        # Plasmid specific exceptions
        if len(complete_assemblies) > 1 and plasmids_only:
            raise AssemblyException('This assembly produces more than one plasmid product. Check your input sequences.')

        if len(complete_assemblies) == 0 and plasmids_only:
            raise AssemblyException('This assembly does not produce any complete products. Check your input sequences.')

        # --- Dump assembly and things into a new Plasmid --- #

        if plasmids_only:

            final_assembly_product = complete_assemblies[0]

            # Find features that actually exist in Plasmid
            plasmid_feature_set = set()
            for feature in final_assembly_product['features']:
                regex_pattern = f'{feature.sequence}|{feature.reverse_complement()}'
                if len(re.findall(regex_pattern, final_assembly_product['sequence'])) > 0:
                    plasmid_feature_set.add(feature)

            return Plasmid(final_assembly_product['sequence'], entity_id=new_id, name=new_id, features=list(plasmid_feature_set),
                           description=new_description, source=final_assembly_product['sources'])


# todo: convert to function
class HomologyAssembly(CloningReaction):
    """
    Perform assemblies that rely on homology complementation methods.

    * Gibson
    * SLIC
    * CPEC/SOE
    """

    def __init__(self, input_dna_list):
        super().__init__(input_dna_list, restriction_enzyme_list=None)

    # --- Setters --- #

    @property
    def restriction_enzyme_list(self):
        return self._restriction_enzyme_list

    @restriction_enzyme_list.setter
    def restriction_enzyme_list(self, input_enzyme_list):
        self._restriction_enzyme_list = input_enzyme_list

    # --- Methods --- #

    def perform_assembly(self, plasmids_only=True, new_id='New_assembly', new_description='New assembly'):
        """
        We're going to use our knowledge of how homology-directed assembly works to take a few short-cuts...
        Nodes are sequences and edges are homology overlaps
        :return:
        """

        homology_overlap = 6

        # --- Assemble Graph --- #

        directed_graph = networkx.MultiDiGraph()

        for part_1, part_2 in combinations(self.input_dna_list, 2):

            print(part_1)
            print(part_2)

            # Get possible part_1 3' overlap with part_2 5'
            part_1_right = part_1.sequence[-homology_overlap:]
            re_pattern = re.compile(f'(?={part_1_right})')
            for match in re_pattern.finditer(part_2.sequence):

                overlap_length = match.start() + len(part_1_right)
                print(part_2.sequence[:overlap_length], part_1.sequence[-overlap_length:])
                if part_2.sequence[:overlap_length] == part_1.sequence[-overlap_length:]:
                    print('Overlap found, adding nodes and edges to graph.')
                    overlap_seq = part_1.sequence[-overlap_length:]

                    directed_graph.add_node(part_1)
                    directed_graph.add_node(part_2)
                    directed_graph.add_edge(part_1, part_2, overlap=overlap_seq)

            # Get possible part_1 5' overlap with part_2 3'
            part_1_left = part_1.sequence[:homology_overlap]
            re_pattern = re.compile(f'(?={part_1_left})')
            for match in re_pattern.finditer(part_2.sequence):
                overlap_length = len(part_2.sequence) - match.start()
                print(part_2.sequence[-overlap_length:], part_1.sequence[:overlap_length])
                if part_2.sequence[-overlap_length:] == part_1.sequence[:overlap_length]:
                    print('Overlap found, adding nodes and edges to graph.')
                    overlap_seq = part_1.sequence[-overlap_length:]

                    directed_graph.add_node(part_1)
                    directed_graph.add_node(part_2)
                    directed_graph.add_edge(part_2, part_1, overlap=overlap_seq)

        # --- Traverse Graph --- #

        # First: simple_cycles
        graph_cycles = networkx.algorithms.cycles.simple_cycles(directed_graph)

        processed_cycles = list()
        all_possible_assemblies = list()

        for cycle in graph_cycles:

            if not any([set(cycle) == processed_cycle for processed_cycle in processed_cycles]):
                sequences_list = list()

                # Iterate through nodes
                for node_l, node_r in pairwise(cycle):
                    current_edge = directed_graph[node_l][node_r]
                    new_sequence_list = [current_edge[index]['overlap'] for index, edge in enumerate(current_edge)]
                    sequences_list.append(new_sequence_list)

                # Get last part of cycle
                last_edge = directed_graph[cycle[-1]][cycle[0]]
                new_sequence_list = [last_edge[index]['overlap'] for index, edge in enumerate(last_edge)]
                sequences_list.append(new_sequence_list)
                all_possible_assemblies += product(*sequences_list)

        # Actually do assemblies
        intermediate_assemblies = list()
        for assembly in all_possible_assemblies:
            print(assembly)
            assembly_dict = dict()
            assembly_dict['sequence'] = ''
            assembly_dict['features'] = list()
            assembly_dict['sources'] = list()
            assembly_dict['description'] = list()

            # Technically it shouldn't be possible for two edges to exist a pair of sequences with homology...
            # but I wanted to reuse sequence_from_cycles()

            for part_1, part_2 in pairwise(assembly):

                overlap_length = len(directed_graph[part_1][part_2][0]['overlap'])
                assembly_dict['sequence'] += part_1.sequence[:-overlap_length]

                if part_1.features:
                    assembly_dict['features'] += part_1.features
                if part_1.description:
                    assembly_dict['description'].append(part_1.description)
                assembly_dict['sources'].append(part_1.source)

            # Get last part
            first_part = assembly[0]
            last_part = assembly[-1]
            overlap_length = len(directed_graph[last_part][first_part][0]['overlap'])
            assembly_dict['sequence'] += last_part.sequence[:-overlap_length]

            if last_part.features:
                assembly_dict['features'] += last_part.features
            if last_part.description:
                assembly_dict['description'].append(last_part.description)
            assembly_dict['sources'].append(last_part.source)

            intermediate_assemblies.append(assembly_dict)

        # --- Make complete assemblies --- #
        """All assemblies are valid since there isn't a digestion step like with sticky-end methods"""
        complete_assemblies = intermediate_assemblies

        # --- Raise exceptions if something went wrong --- #

        # Plasmid specific exceptions
        if len(complete_assemblies) > 1 and plasmids_only:
            raise AssemblyException('This assembly produces more than one plasmid product. Check your input sequences.')

        if len(complete_assemblies) == 0 and plasmids_only:
            raise AssemblyException('This assembly does not produce any complete products. Check your input sequences.')

        # --- Dump assembly and things into a new Plasmid --- #

        if plasmids_only:

            final_assembly_product = complete_assemblies[0]

            # Find features that actually exist in Plasmid
            plasmid_feature_set = set()
            for feature in final_assembly_product['features']:
                regex_pattern = f'{feature.sequence}|{feature.reverse_complement()}'
                if len(re.findall(regex_pattern, final_assembly_product['sequence'])) > 0:
                    plasmid_feature_set.add(feature)

            return Plasmid(final_assembly_product['sequence'], entity_id=new_id, name=new_id, features=list(plasmid_feature_set),
                           description=new_description, source=final_assembly_product['sources'])

class ReactionDefinitionException(Exception):
    pass


class AssemblyException(Exception):
    pass
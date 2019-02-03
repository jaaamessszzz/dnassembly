#! /usr/bin/env python3

import re
import networkx

from .cloning import CloningReaction
from ..dna import Plasmid
from ..utils import pairwise

class GoldenGate(CloningReaction):
    """
    Perfoms golden gate assembly using input parts

    :return:
    """

    def __init__(self, input_dna_list, restriction_enzyme_list):
        super(GoldenGate, self).__init__(input_dna_list, restriction_enzyme_list)

    def verify_parts(self):
        """
        Verify that none of the parts contain an internal restriction site
        :return:
        """
        pass

    def perform_assembly(self, plasmids_only=True, new_id='New_assembly', new_description='New assembly'):
        """
        Performs an assembly provided a pool of Parts.

        Inspiration from pydna to use networkx to represent assemblies as a graph. However I think we can be a lot more
        efficient and enforce our biological constraints at the same time.
        :param plasmids_only: only report cyclic assemblies that fulfill reaction constraints
        :return:
        """
        # --- Assemble Graph --- #

        directed_graph = networkx.DiGraph()

        # Add nodes and edges for each part in digest_pool
        for part in self.digest_pool:

            # Find indicies of sticky ends, continue if both ends to not have sticky ends
            # We will handle this later, this step will mess things up for linear assemblies

            sticky_match_l = part.overhang_5
            sticky_match_r = part.overhang_3

            # Process sticky ends into nodes

            # Nodes store overhang information as ▼TATG▲ or ▲TATG▼
            # I can check {is TATG▲ from Part(TATG▲NNNNNN) in node ▼TATG▲}
            # This way edge are only formed between complementary overhangs

            # (5, TATG) as nodes, where 5 is the overhang strand
            # complementary sticky ends are always formed by the same strand
            # 5'/3' end of DNA is encoded by direction of edges

            if sticky_match_l is None:
                l_node = 'Blunt_Left'
            else:
                l_node = (part.sequence[:sticky_match_l[0]], sticky_match_l[1])

            if sticky_match_r is None:
                r_node = 'Blunt_Right'
            else:
                r_node = (part.sequence[-sticky_match_r[0]:], sticky_match_r[1])  # (overhang_sequence, overhang_strand)

            directed_graph.add_node(l_node)
            directed_graph.add_node(r_node)

            # Add directed edge between nodes (5' -> 3') and add Part as edge attribute
            directed_graph.add_edge(l_node, r_node, {'part': part})

        # --- Traverse Graph --- #

        # First: simple_cycles
        graph_cycles = networkx.algorithms.cycles.simple_cycles(directed_graph)

        # Store already processed cycles as set(nodes)
        processed_cycles = list()
        intermediate_assemblies = list()

        for cycle in graph_cycles:

            if not any([set(cycle) == processed_cycle for processed_cycle in processed_cycles]):
                assembly_dict = dict()

                feature_list = list()
                source_list = list()
                assembled_sequence = ''

                # Iterate through edges and get sequence up to 3' sticky end
                for node_l, node_r in pairwise(cycle):
                    current_part = directed_graph[node_l][node_r]['part']
                    right_sticky_end = 0 if current_part.overhang_3 is None else current_part.overhang_3[0]
                    assembled_sequence = assembled_sequence + current_part.sequence[:-right_sticky_end]
                    feature_list = feature_list + current_part.features
                    source_list.append(current_part.source)

                # Get last part of cycle
                last_part = directed_graph[cycle[-1]][cycle[0]]['part']
                right_sticky_end = 0 if last_part.overhang_3 is None else last_part.overhang_3[0]
                assembled_sequence = assembled_sequence + last_part.sequence[:-right_sticky_end]
                feature_list = feature_list + last_part.features
                source_list.append(last_part.source)

                assembly_dict['sequence'] = assembled_sequence
                assembly_dict['features'] = feature_list
                assembly_dict['sources'] = source_list

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
            plasmid_feature_list = list()
            for feature in final_assembly_product['features']:
                if feature.sequence in final_assembly_product['sequence']:
                    plasmid_feature_list.append(feature)

            return Plasmid(final_assembly_product['sequence'], id=new_id, name=new_id, features=plasmid_feature_list,
                           description=new_description, source=final_assembly_product['sources'])


class AssemblyException(Exception):
    pass
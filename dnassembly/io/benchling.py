# TODO - remove

#! /usr/bin/env python3

from Bio import SeqIO

from ..utils.conversion import convert_dnassembly_to_biopython, convert_biopython_to_dnassembly
from ..utils.benchlingAPI import *

# --- Read files --- #

def read_benchling(benchlingSeqIds, output_format):
    """
    Parse a Benchling Sequence ID and return a DNA entity.

    :param benchlingSeqIds: DNA sequence IDs to pull
    :param output_format: OutputFormat class object for your desired output format
    :return: DNA entity
    """
    parsed_Benchling = getBenchling('dna-sequences/',benchlingSeqId)
    dna_entity = convert_benchling_to_dnassembly(parsed_Benchling, output_format)

    return dna_entity

# --- Write files --- # TODOFIX - could use this inside the Benchling API fn?

'''
def write_benchling(ipnut_dna, output='new_assembly.gb', to_stream=True):
    """
    Output a plasmid in the GenBank file format (.gb)
    :param ipnut_dna: DNA entity
    :param output_path: path
    :return:
    """
    input_as_biopython = convert_dnassembly_to_biopython(ipnut_dna)

    if to_stream:
        SeqIO.write(input_as_biopython, output, "gb")
    else:
        with open(output, "w") as output_path:
            SeqIO.write(input_as_biopython, output_path, "gb")
'''

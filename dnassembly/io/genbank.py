#! /usr/bin/env python3

from Bio import SeqIO

from ..utils.conversion import convert_dnassembly_to_biopython, convert_biopython_to_dnassembly

# --- Read files --- #

def read_genbank(file_path, output_format):
    """
    Parse a GenBank (.gb) file and return a DNA entity. This function assumes that a single GenBank file corresponds to
    one DNA sequence/entity.

    :param file_path: path to .gb file
    :param output_format: OutputFormat class object for your desired output format
    :return: DNA entity
    """
    parsed_genbank = list(SeqIO.parse(file_path, 'gb'))[0]
    dna_entity = convert_biopython_to_dnassembly(parsed_genbank, output_format)

    return dna_entity

# --- Write files --- #

def write_genbank(ipnut_dna, output='new_assembly.gb', to_stream=True):
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
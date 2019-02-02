#! /usr/bin/env python3

from itertools import tee

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation

from ..dna import Feature
from ..utils import pairwise

# --- Read files --- #

def read_genbank(file_path, output_format):
    """
    Parse a GenBank (.gb) file and return a DNA entity. This function assumes that a single GenBank file corresponds to
    one DNA sequence/entity.

    :param file_path: path to .gb file
    :param output_format: OutputFormat class object for your desired output format
    :param input_ape: read unique feature IDs associated with the ApE file format (why would you do this?)
    :return: dict with file contents
    """
    # --- Parse file with BioPython --- #
    parsed_genbank = list(SeqIO.parse(file_path, 'gb'))[0]
    dna_sequence = str(parsed_genbank.seq)

    # todo: automatically infer DNA entity based on annotation
    # I would do it now but all my plasmids are annotated as linear at the moment...

    # Pull any features out of genbank and stick them into Feature class instances
    # Biopython SeqFeature objects only hold feature indicies, which will get lost during an assembly

    feature_list = list()
    feature_sequence_set = set()

    for feature in parsed_genbank.features:

        # Get feature sequence
        if isinstance(feature.location, FeatureLocation):
            feature_sequence = dna_sequence[feature.location.start:feature.location.end]

        # If CompoundLocation, convert feature into regex (e.g. GGTCTC.TATG for a BsaI site + sticky end)
        elif isinstance(feature.location, CompoundLocation):
            feature_parts = feature.location.parts
            indices = sorted([(feature_part.start, feature_part.end) for feature_part in feature_parts])

            # Start with first feature part in feature_sequence
            first_feature_indicies = indices[0]
            feature_sequence = dna_sequence[first_feature_indicies[0]:first_feature_indicies[1]]

            # Add . for nucleotides between features, then add next feature part
            for first_part, second_part in pairwise(indices):
                space_between_features = second_part[0] - first_part[1]
                feature_sequence += '.' * space_between_features
                feature_sequence += dna_sequence[second_part[0]:second_part[1]]

        else:
            print('Only features compatible with BioPython FeatureLocation or CompoundLocation objects are currently supported...')
            continue

        feature_type = feature.type
        strand = feature.location.strand

        # Get feature name... what the difference between locus tag and label is...
        feature_keys = list({'label', 'locus_tag'} & set(feature.qualifiers.keys()))

        if len(feature_keys) == 0:
            feature_name = 'Undefined Feature'
        else:
            feature_name = feature.qualifiers.get(feature_keys[0])[0]

        parsed_feature = Feature(feature_sequence, feature_type, strand, name=feature_name)

        if len(feature_sequence_set & {parsed_feature.reverse_complement(), parsed_feature.sequence}) == 0:
            feature_list.append(parsed_feature)
            feature_sequence_set.add(feature_sequence)

    # Convert BioPython Seq into whatever is defined in output_format
    dna_entity = output_format.value(dna_sequence,
                                     id=parsed_genbank.id,
                                     name=parsed_genbank.name,
                                     description=parsed_genbank.description,
                                     features=feature_list
                                     )

    return dna_entity

# --- Write files --- #

def write_genbank(file_contents, output_path='new_assembly.gb', output_ape=False):
    """
    Output a plasmid in the GenBank file format (.gb)
    :param file_contents: dict containing information to write to file
    :param output_path: path
    :param output_ape: write unique feature IDs associated with the ApE file format (why would you do this?)
    :return:
    """
    pass
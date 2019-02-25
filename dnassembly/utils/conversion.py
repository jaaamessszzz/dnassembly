#! /usr/bin/env python3

import re

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from ..dna import Plasmid, Feature
from ..utils import pairwise, reverse_complement


def convert_biopython_to_dnassembly(parsed_genbank, output_format):
    """
    Convert a BioPython SeqRecord into a DNAssembly DNA entity
    :param seqrecord: BioPython SeqRecord object
    :param output_format: OutputFormat object for your desired output format
    :return:
    """
    dna_sequence = str(parsed_genbank.seq)

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
            print(
                'Only features compatible with BioPython FeatureLocation or CompoundLocation objects are currently supported...')
            continue

        # todo: move all this to alternate Feature constructor that accepts a sequence and BioPython Feature
        feature_type = feature.type
        strand = feature.location.strand

        # get reverse complement of feature sequence if strand == -1
        if strand == -1:
            feature_sequence = reverse_complement(feature_sequence)

        # Get feature name... what the difference between locus tag and label is...
        feature_label_keys = list({'label', 'locus_tag'} & set(feature.qualifiers.keys()))

        if len(feature_label_keys) == 0:
            feature_name = 'Undefined Feature'
        else:
            feature_name = feature.qualifiers.get(feature_label_keys[0])[0]

        # Get feature colors
        forward_color = feature.qualifiers.get('ApEinfo_fwdcolor')[0] if feature.qualifiers.get(
            'ApEinfo_fwdcolor') else None
        reverse_color = feature.qualifiers.get('ApEinfo_revcolor')[0] if feature.qualifiers.get(
            'ApEinfo_revcolor') else None

        # Make Feature
        parsed_feature = Feature(feature_sequence, feature_type, strand, entity_id=feature_name, name=feature_name,
                                 forward_color=forward_color, reverse_color=reverse_color)

        # Add feature to feature_list if feature has not already been encountered
        if len(feature_sequence_set & {parsed_feature.reverse_complement(), parsed_feature.sequence}) == 0:
            feature_list.append(parsed_feature)
            feature_sequence_set.add(feature_sequence)

    # Convert BioPython Seq into whatever is defined in output_format
    dna_entity = output_format.value(dna_sequence,
                                     entity_id=parsed_genbank.id,
                                     name=parsed_genbank.name,
                                     description=parsed_genbank.description,
                                     features=feature_list
                                     )
    return dna_entity


def convert_dnassembly_to_biopython(dna_entity):
    """
    Convert a DNAssembly DNA entity into a BioPython SeqRecord
    :param dna_entity: DNAssembly DNA entity
    :return:
    """
    seq = Seq(dna_entity.sequence, generic_dna)

    # Create feature_map if dna_entity has not been mapped
    if dna_entity.feature_map is None:
        dna_entity.map_features()

    # Convert DNAssembly Features into BioPython Features
    biopython_feature_list = list()
    if dna_entity.feature_map:
        for feature_location in dna_entity.feature_map:
            bound_5, bound_3, strand = feature_location
            for feature in dna_entity.feature_map[feature_location]:
                qualifiers = {"label": feature.name,
                              "ApEinfo_label": feature.name,
                              "ApEinfo_fwdcolor": feature.forward_color,
                              "ApEinfo_revcolor": feature.reverse_color,
                              "ApEinfo_graphicformat": "arrow_data {{0 1 2 0 0 -1} {} 0}"
                              }

                # Process CompoundLocation Features
                if '.' in feature.sequence:
                    compound_feature_list = list()
                    block_start = bound_5

                    # Feature sequence is regex, so find all instances of . to create FeatureLocation blocks
                    compiled_feature = re.compile('\.+')
                    feature_sequence = feature.reverse_complement() if strand == -1 else feature.sequence
                    for match in compiled_feature.finditer(feature_sequence):
                        # Create FeatureLocation block as part of CompoundLocation
                        compound_feature_list.append(FeatureLocation(block_start, block_start + match.start()))
                        block_start = block_start + match.end()
                    # Handle last FeatureLocation block and add to CompoundLocation
                    compound_feature_list.append(FeatureLocation(block_start, bound_3))
                    reverse_order = True if strand == -1 else False
                    compound_feature_list = sorted(compound_feature_list, key=lambda x: x.start, reverse=reverse_order)
                    feature_location = CompoundLocation(compound_feature_list)

                else:
                    feature_location = FeatureLocation(bound_5, bound_3, strand)

                new_seqfeature = SeqFeature(feature_location, strand=strand, id=feature.name, type=feature.feature_type,
                                            qualifiers=qualifiers)
                biopython_feature_list.append(new_seqfeature)

    dna_as_biopython = SeqRecord(seq, id=dna_entity.entity_id, name=dna_entity.name, features=biopython_feature_list, description=dna_entity.description)
    # Add missing annotations
    if isinstance(dna_entity, Plasmid):
        dna_as_biopython.annotations["topology"] = 'circular'
    else:
        dna_as_biopython.annotations["topology"] = 'linear'
    dna_as_biopython.annotations["molecule_type"] = 'DNA'

    return dna_as_biopython
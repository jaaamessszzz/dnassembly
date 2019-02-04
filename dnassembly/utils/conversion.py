#! /usr/bin/env python3

import re

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from ..dna import Plasmid

def convert_biopython_to_dnassembly(seqrecord, output_format):
    """
    Convert a BioPython SeqRecord into a DNAssembly DNA entity
    :param seqrecord: BioPython SeqRecord object
    :param output_format: OutputFormat object for your desired output format
    :return:
    """
    # todo: copy most of genbank.read_genbank code here
    pass

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
                feature_location = CompoundLocation(compound_feature_list)

            else:
                feature_location = FeatureLocation(bound_5, bound_3, strand)

            new_seqfeature = SeqFeature(feature_location, strand=strand, id=feature.name, type=feature.feature_type,
                                        qualifiers=qualifiers)
            biopython_feature_list.append(new_seqfeature)
    # todo: output files annotated with linear/circular depending on dna_entity type... lives somewhere in SeqRecord.feature.qualifiers...

    dna_as_biopython = SeqRecord(seq, id=dna_entity.id, name=dna_entity.name, features=biopython_feature_list, description=dna_entity.description)
    # Add missing annotations
    if isinstance(dna_entity, Plasmid):
        dna_as_biopython.annotations["topology"] = 'circular'
    else:
        dna_as_biopython.annotations["topology"] = 'linear'
    dna_as_biopython.annotations["molecule_type"] = 'DNA'

    return dna_as_biopython
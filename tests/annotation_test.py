from dnassembly.utils.annotation import annotate_moclo
from typing import List, Tuple
import pytest

@pytest.fixture(scope='class')
def moclo_assemblies():
    return [
        ("seq1", "part", []),
        ("seq2", "part", []),
        ("seq3", "cassette", []),
        ("seq4", "cassette", [])
    ]

@pytest.fixture(scope='class')
def invalid_moclo_assemblies():
    return [
        ("seq1", "", []),
        ("seq2", "multicassette", [])
    ]

def test_annotate_moclo(moclo_assemblies: List[Tuple[str, str, List[str]]]):
    """[summary]
    """

    for moclo_assembly in moclo_assemblies:

        dna_sequence = moclo_assembly[0]
        annotation_type = moclo_assembly[1]
        expected_result = moclo_assembly[2]

        actual_result = annotate_moclo(dna_sequence, annotation_type)

        assert len(actual_result) == len(expected_result)
        assert annotate_moclo(dna_sequence, annotation_type) == expected_result

def test_annotate_moclo_exception(invalid_moclo_assemblies: List[Tuple[str, str, List[str]]]):
    """[summary]
    """
    for moclo_assembly in invalid_moclo_assemblies:

        dna_sequence = moclo_assembly[0]
        annotation_type = moclo_assembly[1]
        expected_result = moclo_assembly[2]

        with pytest.raises(Exception):
            actual_result = annotate_moclo(dna_sequence, annotation_type)

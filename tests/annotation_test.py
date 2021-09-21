from dnassembly.utils.annotation import annotate_moclo
from typing import List, Tuple
import pytest

@pytest.fixture(scope='class')
def moclo_assemblies():
    return [
        ("seq1", "part", []),
        ("seq2", "cassette", [])
    ]


def test_annotate_moclo(moclo_assembly: Tuple[str, str, List[str]]):
    """[summary]
    """
    dna_sequence = moclo_assembly[0]
    annotation_type = moclo_assembly[1]
    expected_result = moclo_assembly[2]

    actual_result = annotate_moclo(dna_sequence, annotation_type)

    assert len(actual_result) == len(expected_result)
    assert annotate_moclo(dna_sequence, annotation_type) == expected_result

def test_annotate_moclo_exception():
    """[summary]
    """
    with pytest.raises(Exception):
        annotate_moclo('', 'multicassette')

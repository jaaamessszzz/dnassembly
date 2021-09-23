from dnassembly.utils.annotation import annotate_moclo
from typing import List, Tuple
import pytest

@pytest.fixture(scope='class')
def moclo_assemblies():
    return [
        ("gcatcgtctcattctgaagactcgctagagtagagacgattgaggtcttcaacgcgagacgcgta", "part", ["5", "6"]),
        ("gcatcgtctcattctgaagactcgctagagtagagacgattgaggtcttcaacgcgagacgcgta", "part", ["5", "6"]),
        #("", "part", []), TODO - configure function to return an empty list instead of a bool
        #("seq3", "cassette", []),
        #("seq4", "cassette", [])
    ]

def test_annotate_moclo(moclo_assemblies: List[Tuple[str, str, List[str]]]):
    """
    Given a sequence, annotation type, and expected annotations 
    assert that the actual annotation equals the expected annotation
    """

    for moclo_assembly in moclo_assemblies:

        dna_sequence = moclo_assembly[0]
        annotation_type = moclo_assembly[1]
        expected_result = moclo_assembly[2]

        actual_result = annotate_moclo(dna_sequence, annotation_type)

        actual_result.sort()
        expected_result.sort()

        # assert that actual and expected results output the same number of parts
        assert len(actual_result) == len(expected_result)
        # assert that both lists have the same parts
        for idx in range(len(actual_result)):
            assert actual_result[idx] == expected_result[idx]

@pytest.fixture(scope='class')
def invalid_moclo_assemblies():
    return [
        ("atat", ""),
        ("atat", "multicassette"),
        ("", "")
    ]

def test_annotate_moclo_exception(invalid_moclo_assemblies: List[Tuple[str, str]]):
    """
    Asserts that annotate_moclo() will throw an error with invalid input parameters
    """
    for moclo_assembly in invalid_moclo_assemblies:

        dna_sequence = moclo_assembly[0]
        annotation_type = moclo_assembly[1]

        with pytest.raises(Exception):
            annotate_moclo(dna_sequence, annotation_type)

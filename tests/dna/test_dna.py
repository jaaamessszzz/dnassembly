from dnassembly.dna import DNA, SequenceException, RestrictionSiteException
from typing import List, Tuple
import pytest

#TODO - create test_dna_constructor
#TODO - create test_dna_repr
#TODO - create test_dna_strand_setter
#TODO - create test_dna_strand
#TODO - create test_dna_map_features
#TODO - create test_dna_find_cut_indicies

@pytest.fixture(scope='class')
def valid_dna_sequences() -> List[Tuple[str, str]]:
    return [
        ("ATAC", "TATG"),
        ("ATAC", "TATG"),
        ("AAA", "TTT"),
        ("TTT", "AAA"),
        ("ATcGACaTtAcc", "TAGCTGTAATGG"),
        ("ATGGACATTAGG", "TACCTGTAATCC"),
        ("GGACCATGAC", "CCTGGTACTG"),
        ("GGaCcATcca", "CCTGGTAGGT")
    ]

@pytest.fixture(scope='class')
def invalid_dna_sequences() -> List[str]:
    return [
        "",
        "nnn",
        "aaan",
        "TTT+",
        "ATGGACATTAGGN",
        "gGACCaTGacn",
        "132"
    ]

def test_dna_seq_setter(valid_dna_sequences: List[Tuple[str, str]]):
    """
    Tests the setter method for valid sequences
    """
    dna = DNA("AAA")
    
    for dna_sequence_data in valid_dna_sequences:

        sequence = dna_sequence_data[0]
        dna.sequence = sequence
        assert dna.sequence == sequence.upper()

def test_dna_seq_setter_invalid(invalid_dna_sequences: List[str]):
    """
    Tests the setter method for invalid sequences
    """
    dna = DNA("AAA")

    for sequence in invalid_dna_sequences:
        
        with pytest.raises(SequenceException):
            print(sequence)
            dna.sequence = sequence

def test_dna_len(valid_dna_sequences: List[Tuple[str, str]]):
    """
    Tests the length class method for DNA
    """
    for dna_sequence_data in valid_dna_sequences:

        sequence = dna_sequence_data[0]
        expected_seq_length = len(sequence)
        assert len(DNA(sequence)) == expected_seq_length

def test_dna_complement(valid_dna_sequences: List[Tuple[str, str]]):
    """
    Tests the complement function for DNA
    """
    for dna_sequence_data in valid_dna_sequences:

        sequence = dna_sequence_data[0]
        comp_sequence = dna_sequence_data[1]
        assert DNA(sequence).complement() == comp_sequence

def test_dna_reverse_complement(valid_dna_sequences: List[Tuple[str, str]]):
    """
    Tests the reverse_complement function for DNA
    """
    for dna_sequence_data in valid_dna_sequences:

        sequence = dna_sequence_data[0]
        comp_sequence = dna_sequence_data[1]
        reverse_comp_sequence = comp_sequence[::-1]
        assert DNA(sequence).reverse_complement() == reverse_comp_sequence

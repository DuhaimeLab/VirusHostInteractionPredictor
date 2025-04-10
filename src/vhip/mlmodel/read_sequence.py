"""Utility functions to read fasta sequences.

This module provides:
- read_sequence: read f asta file using path as input
- read_headers: retrieve the headers for fasta path
"""


import json
from typing import Any, Dict, List

from Bio import SeqIO  # pyright: ignore[reportMissingTypeStubs]

complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
}


def read_sequence(path: str) -> list[str]:
    """Return the sequence in a fasta file as a string.

    Args:
        path (str): path to fasta file

    Returns:
        str: nucleotide sequence
    """
    sequences: List[str] = []
    for record in SeqIO.parse(path, "fasta"):  # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
        sequences.append(str(record.seq))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    return sequences


def read_headers(path: str) -> list[str]:
    """Return all headers in fasta file as a list.

    Args:
        path (str): path to fasta file

    Returns:
        list: list of all headers for a given fasta file
    """
    result = []
    for record in SeqIO.parse(path, "fasta"):  # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
        result.append(record.id)  # pyright: ignore[reportUnknownMemberType]
    return result  # pyright: ignore[reportUnknownVariableType]


def read_annotated_genes(path: str) -> List[Dict[str, Any]]:
    """Return all components of a .json fasta file output from bakta.

    Args:
        path (str): path to annotated gene json file

    Returns:
        List[Dict[str,Any]]: list of all features (i.e. genes) in the input .json file
    """
    with open(path, 'r') as file:
        annotations = json.load(file)

    features: List[Dict[str, Any]] = annotations['features']
    return features


def reverse_complement(sequence: str) -> str:
    """Get reverse complement of a sequence (for work with anticodons).

    Args:
        sequence (str): nucleotide sequence of any length, with upper or lower case letters

    Returns:
        reverse_complement_sequence: reverse complement of the input sequence, in upper case letters
    """
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = "".join(
        complement.get(base, base) for base in reverse_sequence
    )
    return reverse_complement_sequence.upper()

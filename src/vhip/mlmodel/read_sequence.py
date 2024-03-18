"""Utility functions to read fasta sequences.

This module provides:
- read_sequence: read f asta file using path as input
- read_headers: retrieve the headers for fasta path
"""

from typing import List

from Bio import SeqIO  # pyright: ignore[reportMissingTypeStubs]


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

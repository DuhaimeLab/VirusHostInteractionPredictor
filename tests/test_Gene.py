"""Pytest for Gene methods in gene_features module."""

from vhip.mlmodel.gene_features import Gene


def test_Gene_attributes():
    """Test code to create Gene object and initialize class attributes."""
    test_gene = Gene("ATGCCGATT")
    assert test_gene.seq == "ATGCCGATT"
    assert test_gene.codon_length == 3



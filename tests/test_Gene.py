"""Pytest for Gene methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import Gene


def test_Gene_init():
    """Test code to create Gene object and initialize class attributes."""
    test_gene_1 = Gene(gene_seq="ATGCCGATT",gene_id="test_gene_1",gene_product="test_gene_product_1")
    assert test_gene_1.seq == "ATGCCGATT"
    assert test_gene_1.codon_length == 3
    assert test_gene_1.gene_id == "test_gene_1"
    assert test_gene_1.gene_product == "test_gene_product_1"

    # Test that Exception is raised when gene length not divisible by codon length
    with pytest.raises(Exception):
        Gene("ATGCCGATTA")


def test_Gene_calculate_codon_counts():
    """Test code to calculate codon counts for a given gene."""
    test_gene = Gene("NTGCCGATT")
    test_gene.calculate_codon_counts()
    assert test_gene.codon_dict == {
        "ATA": 0,
        "ATC": 0,
        "ATT": 1,
        "ATG": 0,
        "ACA": 0,
        "ACC": 0,
        "ACG": 0,
        "ACT": 0,
        "AAC": 0,
        "AAT": 0,
        "AAA": 0,
        "AAG": 0,
        "AGC": 0,
        "AGT": 0,
        "AGA": 0,
        "AGG": 0,
        "CTA": 0,
        "CTC": 0,
        "CTG": 0,
        "CTT": 0,
        "CCA": 0,
        "CCC": 0,
        "CCG": 1,
        "CCT": 0,
        "CAC": 0,
        "CAT": 0,
        "CAA": 0,
        "CAG": 0,
        "CGA": 0,
        "CGC": 0,
        "CGG": 0,
        "CGT": 0,
        "GTA": 0,
        "GTC": 0,
        "GTG": 0,
        "GTT": 0,
        "GCA": 0,
        "GCC": 0,
        "GCG": 0,
        "GCT": 0,
        "GAC": 0,
        "GAT": 0,
        "GAA": 0,
        "GAG": 0,
        "GGA": 0,
        "GGC": 0,
        "GGG": 0,
        "GGT": 0,
        "TCA": 0,
        "TCC": 0,
        "TCG": 0,
        "TCT": 0,
        "TTC": 0,
        "TTT": 0,
        "TTA": 0,
        "TTG": 0,
        "TAC": 0,
        "TAT": 0,
        "TAA": 0,
        "TAG": 0,
        "TGC": 0,
        "TGT": 0,
        "TGA": 0,
        "TGG": 0,
    }
    assert test_gene.number_imprecise_codons == 1


def test_Gene_calculate_aa_counts():
    """Test code to calculate amino acid counts for a given gene."""
    test_gene = Gene("NTGCCGATT")
    test_gene.calculate_aa_counts()
    assert test_gene.aa_dict == {
        "_": 0,
        "C": 0,
        "Y": 0,
        "W": 0,
        "S": 0,
        "R": 0,
        "N": 0,
        "F": 0,
        "V": 0,
        "M": 0,
        "I": 1,
        "L": 0,
        "E": 0,
        "A": 0,
        "H": 0,
        "Q": 0,
        "K": 0,
        "D": 0,
        "G": 0,
        "T": 0,
        "P": 1,
    }

"""Pytest for CodonBiasComparison methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import CodonBiasComparison


def test_CodonBiasComparison_init():
    """Test code to create CodonBiasComparison object and initialize class attributes."""
    # Create host GeneSet object and calculate codon bias
    test_host_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_host_GeneSet.codon_counts()
    test_host_GeneSet.codon_frequency()
    test_host_GeneSet.amino_acid_counts()
    test_host_GeneSet.amino_acid_frequency()
    test_host_GeneSet.RSCU()

    # Create virus GeneSet object and calculate codon bias
    test_virus_GeneSet = GeneSet("tests/datatests/test_virus_short_genes_file.ffn")
    test_virus_GeneSet.codon_counts()
    test_virus_GeneSet.codon_frequency()
    test_virus_GeneSet.amino_acid_counts()
    test_virus_GeneSet.amino_acid_frequency()
    test_virus_GeneSet.RSCU()


"""Pytest for CodonBiasComparison methods in gene_features module."""

import scipy

from vhip.mlmodel.gene_features import CodonBiasComparison, GeneSet


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

    # test 1 - test CodonBiasComparsion initializes with input codon counts (GeneSet.codon_dict)
    test_comparison_1 = CodonBiasComparison(test_host_GeneSet.codon_dict, test_virus_GeneSet.codon_dict)
    assert len(test_comparison_1.host_list) == 64
    assert len(test_comparison_1.virus_list) == 64

    # test 2 - test CodonBiasComparsion initializes with input codon frequency (GeneSet.codon_frq)
    test_comparison_2 = CodonBiasComparison(test_host_GeneSet.codon_frq, test_virus_GeneSet.codon_frq)
    assert len(test_comparison_2.host_list) == 64
    assert len(test_comparison_2.virus_list) == 64

    # test 3 - test CodonBiasComparsion initializes with input amino acid counts (GeneSet.aa_dict)
    test_comparison_3 = CodonBiasComparison(test_host_GeneSet.aa_dict, test_virus_GeneSet.aa_dict)
    assert len(test_comparison_3.host_list) == 20
    assert len(test_comparison_3.virus_list) == 20

    # test 4 - test CodonBiasComparsion initializes with input amino acid frequency (GeneSet.aa_frq)
    test_comparison_4 = CodonBiasComparison(test_host_GeneSet.aa_frq, test_virus_GeneSet.aa_frq)
    assert len(test_comparison_4.host_list) == 20
    assert len(test_comparison_4.virus_list) == 20

    # test 5 - test CodonBiasComparsion initializes with input RSCU (GeneSet.RSCU_dict)
    test_comparison_5 = CodonBiasComparison(test_host_GeneSet.RSCU_dict, test_virus_GeneSet.RSCU_dict)
    assert len(test_comparison_5.host_list) == 64
    assert len(test_comparison_5.virus_list) == 64

def test_CodonBiasComparison_slope():
    """Test code to calculate slope between codon biases of a virus GeneSet and host GeneSet."""
    # test 1 - test that slope = 1 when virus and host inputs are exactly the same
    test_host_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_host_GeneSet.codon_counts()
    test_virus_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_virus_GeneSet.codon_counts()
    test_comparison = CodonBiasComparison(test_host_GeneSet.codon_dict, test_virus_GeneSet.codon_dict)
    assert test_comparison.slope() == 1

    # test 2 - test slope calculation when virus and host inputs have differences
    test_host_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_host_GeneSet.codon_counts()
    test_virus_GeneSet = GeneSet("tests/datatests/test_virus_short_genes_file.ffn")
    test_virus_GeneSet.codon_counts()
    test_comparison = CodonBiasComparison(test_host_GeneSet.codon_dict, test_virus_GeneSet.codon_dict)
    expected_host_list = [1] * 4 + [0] * 60
    expected_virus_list = [2] * 2 + [1] * 3 + [0] * 59
    assert test_comparison.slope() == scipy.stats.linregress(expected_host_list, expected_virus_list)[0]



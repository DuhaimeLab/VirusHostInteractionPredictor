"""Pytest for CodonBiasComparison methods in gene_features module."""

import scipy  # type: ignore

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
    assert len(test_comparison_1.host_list) == 64 # there are 64 possible unique codons
    assert len(test_comparison_1.virus_list) == 64

    # test 2 - test CodonBiasComparsion initializes with input codon frequency (GeneSet.codon_frq)
    test_comparison_2 = CodonBiasComparison(test_host_GeneSet.codon_frq, test_virus_GeneSet.codon_frq)
    assert len(test_comparison_2.host_list) == 64
    assert len(test_comparison_2.virus_list) == 64

    # test 3 - test CodonBiasComparsion initializes with input amino acid counts (GeneSet.aa_dict)
    test_comparison_3 = CodonBiasComparison(test_host_GeneSet.aa_dict, test_virus_GeneSet.aa_dict)
    assert len(test_comparison_3.host_list) == 21 # there are 20 possible unique amino acids plus a stop signal
    assert len(test_comparison_3.virus_list) == 21

    # test 4 - test CodonBiasComparsion initializes with input amino acid frequency (GeneSet.aa_frq)
    test_comparison_4 = CodonBiasComparison(test_host_GeneSet.aa_frq, test_virus_GeneSet.aa_frq)
    assert len(test_comparison_4.host_list) == 21
    assert len(test_comparison_4.virus_list) == 21

    # test 5 - test CodonBiasComparsion initializes with input RSCU (GeneSet.RSCU_dict)
    test_comparison_5 = CodonBiasComparison(test_host_GeneSet.RSCU_dict, test_virus_GeneSet.RSCU_dict)
    assert len(test_comparison_5.host_list) == 64
    assert len(test_comparison_5.virus_list) == 64


def test_CodonBiasComparison_methods():
    """Test code to calculate slope, R2, and cosine_similarity between codon biases of a virus GeneSet and host GeneSet.

    Here we use GeneSet.codon_dict (codon counts) dictionaries as the test inputs.
    """
    # test 1 - test that metrics all = 1 when virus and host inputs are exactly the same
    test_host_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_host_GeneSet.codon_counts()
    test_virus_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_virus_GeneSet.codon_counts()

    test_comparison = CodonBiasComparison(test_host_GeneSet.codon_dict, test_virus_GeneSet.codon_dict)
    test_comparison.slope()
    test_comparison.R2()
    test_comparison.cosine_similarity()

    assert test_comparison.slope == 1
    assert test_comparison.R2 == 1
    assert test_comparison.cosine_similarity == 1


    # test 2 - test metrics calculation when virus and host inputs have differences
    test_host_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")
    test_host_GeneSet.codon_counts()
    test_virus_GeneSet = GeneSet("tests/datatests/test_virus_short_genes_file.ffn")
    test_virus_GeneSet.codon_counts()

    test_comparison = CodonBiasComparison(test_host_GeneSet.codon_dict, test_virus_GeneSet.codon_dict)
    test_comparison.slope()
    test_comparison.R2()
    test_comparison.cosine_similarity()

    expected_virus_codon_dict = {
            "ATA": 0,
            "ATC": 0,
            "ATT": 0,
            "ATG": 2,
            "ACA": 0,
            "ACC": 0,
            "ACG": 0,
            "ACT": 0,
            "AAC": 0,
            "AAT": 0,
            "AAA": 1,
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
            "CCG": 0,
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
            "GAA": 2,
            "GAG": 0,
            "GGA": 0,
            "GGC": 0,
            "GGG": 0,
            "GGT": 0,
            "TCA": 1,
            "TCC": 1,
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
            "TGG": 1,
        }
    expected_virus_list = list(expected_virus_codon_dict.values())
    expected_host_codon_dict = {
            "ATA": 0,
            "ATC": 0,
            "ATT": 0,
            "ATG": 1,
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
            "CCG": 0,
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
            "GAA": 1,
            "GAG": 0,
            "GGA": 0,
            "GGC": 0,
            "GGG": 0,
            "GGT": 0,
            "TCA": 1,
            "TCC": 1,
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
    expected_host_list = list(expected_host_codon_dict.values())


    assert test_host_GeneSet.codon_dict == expected_host_codon_dict
    assert test_virus_GeneSet.codon_dict == expected_virus_codon_dict

    assert test_comparison.slope == scipy.stats.linregress(expected_host_list, expected_virus_list)[0]
    assert test_comparison.R2 == scipy.stats.linregress(expected_host_list, expected_virus_list)[2] ** 2
    assert test_comparison.cosine_similarity == 1 - scipy.spatial.distance.cosine(expected_host_list, expected_virus_list)

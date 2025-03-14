"""Pytest for GeneSet methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import GeneSet


def test_GeneSet_init():
    """Test code to create GeneSet object and initialize class attributes."""
    # test 1 - test GeneSet object creation generates expected number of genes in genes attribute
    test_GeneSet = GeneSet("tests/datatests/test_annotated_genes.ffn")
    assert len(test_GeneSet.genes) == 3

    # test 2 - test GeneSet object creation generates expected gene attributes for a gene
    assert (
        test_GeneSet.genes[2].seq
        == "TTGGTTGAAGAAGTAGTTGTAGATGGCGACATCACATTAGGACAATTTCTAAAGACGGAAGGTATTATCGAATCTGGCGGGCAAGCGAAATGGTTCTTAAATGAGTTTGAAGTATTGTTAAACAATACGCGTGAAACACGCCGTGGTAAAAAGTTAAGCCATCGTGACACAATTGAGATACCAGAAATACCTGAAGTGGGTTCATTTGTGATTTTGCATCAAGGTGAAGAATGA"
    )
    assert test_GeneSet.genes[2].gene_id == "ABDEAL_00015"
    assert test_GeneSet.genes[2].gene_product == "S4 domain-containing protein YaaA"

    # test 3 - test Exception is raised if gene file is empty
    with pytest.raises(Exception):
        GeneSet("tests/datatests/test_empty_file.ffn")

    # test 4 - test individual genes are skipped and exception is handled when gene length not divisible by codon length
    test_GeneSet2 = GeneSet("tests/datatests/test_short_genes_file.ffn")
    assert len(test_GeneSet2.genes) == 2
    assert test_GeneSet2.skipped_genes == ["ABDEAL_00010"]


def test_GeneSet_codon_counts():
    """Test code to calculate codon counts across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")

    # test 1 - high threshold_imprecise and high threshold_skipped_genes
    test_GeneSet.codon_counts(1000, 1000)
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 0
    assert test_GeneSet.codon_dict == {
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
        "TCA": 2,
        "TCC": 1,
        "TCG": 0,
        "TCT": 0,
        "TTC": 0,
        "TTT": 0,
        "TTA": 0,
        "TTG": 1,
        "TAC": 0,
        "TAT": 0,
        "TAA": 0,
        "TAG": 0,
        "TGC": 0,
        "TGT": 0,
        "TGA": 0,
        "TGG": 0,
    }

    # test 2 - low threshold_imprecise and high threshold_skipped_genes
    test_GeneSet.codon_counts(0, 1000)
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.codon_dict == {
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

    # test 3 - low threshold threshold_skipped_genes, raising error
    with pytest.raises(Exception):
        test_GeneSet.codon_counts(0, 0)


def test_GeneSet_codon_frequency():
    """Test code to calculate codon frequency across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")

    test_GeneSet.codon_frequency()  # using default thresholds for tolerable imprecise codons and skipped genes
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_genes) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.codon_frq == {
        "ATA": 0.0,
        "ATC": 0.0,
        "ATT": 0.0,
        "ATG": 1 / 4,
        "ACA": 0.0,
        "ACC": 0.0,
        "ACG": 0.0,
        "ACT": 0.0,
        "AAC": 0.0,
        "AAT": 0.0,
        "AAA": 0.0,
        "AAG": 0.0,
        "AGC": 0.0,
        "AGT": 0.0,
        "AGA": 0.0,
        "AGG": 0.0,
        "CTA": 0.0,
        "CTC": 0.0,
        "CTG": 0.0,
        "CTT": 0.0,
        "CCA": 0.0,
        "CCC": 0.0,
        "CCG": 0.0,
        "CCT": 0.0,
        "CAC": 0.0,
        "CAT": 0.0,
        "CAA": 0.0,
        "CAG": 0.0,
        "CGA": 0.0,
        "CGC": 0.0,
        "CGG": 0.0,
        "CGT": 0.0,
        "GTA": 0.0,
        "GTC": 0.0,
        "GTG": 0.0,
        "GTT": 0.0,
        "GCA": 0.0,
        "GCC": 0.0,
        "GCG": 0.0,
        "GCT": 0.0,
        "GAC": 0.0,
        "GAT": 0.0,
        "GAA": 1 / 4,
        "GAG": 0.0,
        "GGA": 0.0,
        "GGC": 0.0,
        "GGG": 0.0,
        "GGT": 0.0,
        "TCA": 1 / 4,
        "TCC": 1 / 4,
        "TCG": 0.0,
        "TCT": 0.0,
        "TTC": 0.0,
        "TTT": 0.0,
        "TTA": 0.0,
        "TTG": 0.0,
        "TAC": 0.0,
        "TAT": 0.0,
        "TAA": 0.0,
        "TAG": 0.0,
        "TGC": 0.0,
        "TGT": 0.0,
        "TGA": 0.0,
        "TGG": 0.0,
    }


def test_GeneSet_amino_acid_counts():
    """Test code to calculate amino acid counts across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")

    test_GeneSet.amino_acid_counts()  # using default thresholds for tolerable imprecise codons and skipped genes
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_genes) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.aa_dict == {
        "I": 0,
        "M": 1,
        "T": 0,
        "N": 0,
        "K": 0,
        "S": 2,
        "R": 0,
        "L": 0,
        "P": 0,
        "H": 0,
        "Q": 0,
        "V": 0,
        "A": 0,
        "D": 0,
        "E": 1,
        "G": 0,
        "F": 0,
        "Y": 0,
        "C": 0,
        "W": 0,
    }


def test_GeneSet_amino_acid_frequency():
    """Test code to calculate amino acid frequency across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")

    test_GeneSet.amino_acid_frequency()  # using default thresholds for tolerable imprecise codons and skipped genes
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_genes) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.aa_frq == {
        "I": 0.0,
        "M": 1 / 4,
        "T": 0.0,
        "N": 0.0,
        "K": 0.0,
        "S": 2 / 4,
        "R": 0.0,
        "L": 0.0,
        "P": 0.0,
        "H": 0.0,
        "Q": 0.0,
        "V": 0.0,
        "A": 0.0,
        "D": 0.0,
        "E": 1 / 4,
        "G": 0.0,
        "F": 0.0,
        "Y": 0.0,
        "C": 0.0,
        "W": 0.0,
    }


def test_GeneSet_RSCU():
    """Test code to calculate relative synonymous codon usage (RSCU) across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes_file.ffn")

    test_GeneSet.RSCU()  # using default thresholds for tolerable imprecise codons and skipped genes
    assert test_GeneSet.imprecise_codons == 1
    assert len(test_GeneSet.skipped_genes) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.RSCU_dict == {
        "ATA": 0.0,
        "ATC": 0.0,
        "ATT": 0.0,
        "ATG": 1 / (1 / 1),
        "ACA": 0.0,
        "ACC": 0.0,
        "ACG": 0.0,
        "ACT": 0.0,
        "AAC": 0.0,
        "AAT": 0.0,
        "AAA": 0.0,
        "AAG": 0.0,
        "AGC": 0.0,
        "AGT": 0.0,
        "AGA": 0.0,
        "AGG": 0.0,
        "CTA": 0.0,
        "CTC": 0.0,
        "CTG": 0.0,
        "CTT": 0.0,
        "CCA": 0.0,
        "CCC": 0.0,
        "CCG": 0.0,
        "CCT": 0.0,
        "CAC": 0.0,
        "CAT": 0.0,
        "CAA": 0.0,
        "CAG": 0.0,
        "CGA": 0.0,
        "CGC": 0.0,
        "CGG": 0.0,
        "CGT": 0.0,
        "GTA": 0.0,
        "GTC": 0.0,
        "GTG": 0.0,
        "GTT": 0.0,
        "GCA": 0.0,
        "GCC": 0.0,
        "GCG": 0.0,
        "GCT": 0.0,
        "GAC": 0.0,
        "GAT": 0.0,
        "GAA": 1 / (1 / 2),
        "GAG": 0.0,
        "GGA": 0.0,
        "GGC": 0.0,
        "GGG": 0.0,
        "GGT": 0.0,
        "TCA": 1 / (2 / 6),
        "TCC": 1 / (2 / 6),
        "TCG": 0.0,
        "TCT": 0.0,
        "TTC": 0.0,
        "TTT": 0.0,
        "TTA": 0.0,
        "TTG": 0.0,
        "TAC": 0.0,
        "TAT": 0.0,
        "TAA": 0.0,
        "TAG": 0.0,
        "TGC": 0.0,
        "TGT": 0.0,
        "TGA": 0.0,
        "TGG": 0.0,
    }


def test_GeneSet_tRNA_counts():
    """Test code to calculate tRNA gene copy counts across a GeneSet."""
    test_GeneSet = GeneSet("tests/datatests/test_host_tRNA_genes.ffn")
    test_GeneSet.tRNA_counts()
    assert test_GeneSet.total_tRNA == 6
    assert test_GeneSet.tRNA_dict_aa == {
        "M": 1,
        "W": 0,
        "I": 0,
        "T": 0,
        "N": 0,
        "K": 0,
        "S": 3,
        "R": 2,
        "L": 0,
        "P": 0,
        "H": 0,
        "Q": 0,
        "V": 0,
        "A": 0,
        "D": 0,
        "E": 0,
        "G": 0,
        "F": 0,
        "Y": 0,
        "C": 0,
    }
    assert test_GeneSet.tRNA_dict_tcc == {
        "ATG": 1,
        "TGG": 0,
        "ATA": 0,
        "ATC": 0,
        "ATT": 0,
        "ACA": 0,
        "ACC": 0,
        "ACG": 0,
        "ACT": 0,
        "AAC": 0,
        "AAT": 0,
        "AAA": 0,
        "AAG": 0,
        "AGC": 2,
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
        "CGT": 2,
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
        "TCC": 1,
        "TCG": 0,
        "TCT": 0,
        "TTC": 0,
        "TTT": 0,
        "TTA": 0,
        "TTG": 0,
        "TAC": 0,
        "TAT": 0,
        "TGC": 0,
        "TGT": 0,
    }


def test_GeneSet_tRNA_frequency():
    """Test code to calculate tRNA gene copy frequency out of total tRNA genes in a GeneSet."""
    test_GeneSet = GeneSet("tests/datatests/test_host_tRNA_genes.ffn")
    test_GeneSet.tRNA_frequency()
    assert test_GeneSet.total_tRNA == 6
    assert test_GeneSet.tRNA_frq_aa == {
        "M": 1 / 6,
        "W": 0,
        "I": 0,
        "T": 0,
        "N": 0,
        "K": 0,
        "S": 3 / 6,
        "R": 2 / 6,
        "L": 0,
        "P": 0,
        "H": 0,
        "Q": 0,
        "V": 0,
        "A": 0,
        "D": 0,
        "E": 0,
        "G": 0,
        "F": 0,
        "Y": 0,
        "C": 0,
    }
    assert test_GeneSet.tRNA_frq_tcc == {
        "ATG": 1 / 6,  # Methionine
        "TGG": 0,
        "ATA": 0,
        "ATC": 0,
        "ATT": 0,
        "ACA": 0,
        "ACC": 0,
        "ACG": 0,
        "ACT": 0,
        "AAC": 0,
        "AAT": 0,
        "AAA": 0,
        "AAG": 0,
        "AGC": 2 / 6,  # Serine
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
        "CGT": 2 / 6,  # Arginine
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
        "TCC": 1 / 6,  # Serine
        "TCG": 0,
        "TCT": 0,
        "TTC": 0,
        "TTT": 0,
        "TTA": 0,
        "TTG": 0,
        "TAC": 0,
        "TAT": 0,
        "TGC": 0,
        "TGT": 0,
    }

"""Pytest for Gene methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import Gene, CDSGene, tRNAGene


def test_Gene_init():
    """Test code to create Gene object and initialize class attributes."""
    # Test 1: Gene with valid sequence and inputs
    test_gene_1 = Gene({"type": "cds", "id": "1", "gene": "test_gene_1", "product": "test_gene_product_1", "nt": "ATGCCGATTTAG", "aa": "MPI"})
    assert test_gene_1.input_error is False
    assert test_gene_1.type == "cds"
    assert test_gene_1.id == "1"
    assert test_gene_1.gene == "test_gene_1"
    assert test_gene_1.product == "test_gene_product_1"

    # Test 2: Gene with generic missing arguments
    test_gene_2 = Gene({"random": "random"})
    assert test_gene_2.input_error is True


def test_CDSGene_init():
    """Test code to create CDSGene object and initialize class attributes."""
    # Test 1: CDS gene with all expected inputs
    test_cds_1 = CDSGene(json_dict={"type": "cds", "id": "1", "gene": "test_gene_1", "product": "test_gene_product_1", "nt": "ATGCCGATTTAG", "aa": "MPI"})
    assert test_cds_1.input_error is False
    assert test_cds_1.cds_len_error is False
    assert test_cds_1.aa_input_error is False
    assert test_cds_1.type == "cds"
    assert test_cds_1.id == "1"
    assert test_cds_1.gene == "test_gene_1"
    assert test_cds_1.product == "test_gene_product_1"
    assert test_cds_1.nt == "ATGCCGATTTAG"
    assert test_cds_1.aa == "MPI"
    assert test_cds_1.codon_length == 3

    # Test 2: Gene with generic missing arguments
    test_gene_2 = CDSGene({"random": "random"})
    assert test_gene_2.input_error is True

    # Test 3: Type not 'cds'
    with pytest.raises(Exception):
        CDSGene({"type": "tRNA", "id": "3", "gene": "test_gene_3", "product": "test_gene_product_3", "nt": "ATGCCGATTTAG", "aa": "MPI"})

    # Test 4: CDS gene with invalid sequence (not divisible by 3)
    test_gene_4 = CDSGene({"type": "cds", "id": "4", "gene": "test_gene_4", "product": "test_gene_product_4", "nt": "ATGCCGATTTAGG", "aa": "MPI"})
    assert test_gene_4.input_error is False
    assert test_gene_4.aa_input_error is False
    assert test_gene_4.cds_len_error is True

    # Test 5: CDS gene with missing aa argument
    test_gene_5 = CDSGene({"type": "cds", "id": "4", "gene": "test_gene_4", "product": "test_gene_product_4", "nt": "ATGCCGATTTAG"})
    assert test_gene_5.aa_input_error is True

def test_CDSGene_calculate_codon_counts():
    """Test code to calculate codon counts for a given CDS gene."""
    # Test 1: CDS gene with 1 imprecise codon
    test_gene = CDSGene({"type": "cds", "id": "1", "gene": "test_gene_1", "product": "test_gene_product_1", "nt": "NTGCCGATTTAG", "aa": "PI"})
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
        "TAG": 1,
        "TGC": 0,
        "TGT": 0,
        "TGA": 0,
        "sTGA": 0,
        "TGG": 0,
    }
    assert len(test_gene.imprecise_codons) == 1

    # Test 2: CDS gene with no imprecise codons
    test_gene_2 = CDSGene({"type": "cds", "id": "2", "gene": "test_gene_2", "product": "test_gene_product_2", "nt": "ATGCCGATTTAG", "aa": "MPI"})
    test_gene_2.calculate_codon_counts()
    assert len(test_gene_2.imprecise_codons) == 0

    # Test 3: CDS gene with a selenocysteine-encoding codon and no TGA stop codon
    test_gene3 = CDSGene({"type": "cds", "id": "3", "gene": "test_gene_3", "product": "test_gene_product_3", "nt": "ATGTGATAA", "aa": "MU"})
    test_gene3.calculate_codon_counts()
    assert test_gene3.codon_dict["TGA"] == 0
    assert test_gene3.codon_dict["sTGA"] == 1

    # Test 4: CDS gene with a selenocysteine-encoding codon and a TGA stop codon
    test_gene4 = CDSGene({"type": "cds", "id": "4", "gene": "test_gene_4", "product": "test_gene_product_4", "nt": "ATGTGATGA", "aa": "MU"})
    test_gene4.calculate_codon_counts()
    assert test_gene4.codon_dict["TGA"] == 1
    assert test_gene4.codon_dict["sTGA"] == 1

    # Test 5: CDS gene with a TGA stop codon but no selenocysteine-encoding codon
    test_gene5 = CDSGene({"type": "cds", "id": "5", "gene": "test_gene_5", "product": "test_gene_product_5", "nt": "ATGTGA", "aa": "M"})
    test_gene5.calculate_codon_counts()
    assert test_gene5.codon_dict["TGA"] == 1
    assert test_gene5.codon_dict["sTGA"] == 0


def test_CDSGene_calculate_aa_counts():
    """Test code to calculate amino acid counts for a given gene."""
    # Test 1: CDS gene with expected amino acids in input dict
    test_gene = CDSGene({"type": "cds", "id": "1", "gene": "test_gene_1", "product": "test_gene_product_1", "nt": "ATGCCGATTTAG", "aa": "MPI"})
    test_gene.calculate_aa_counts()
    assert test_gene.aa_dict == {
        "C": 0,
        "Y": 0,
        "W": 0,
        "S": 0,
        "R": 0,
        "N": 0,
        "F": 0,
        "V": 0,
        "M": 1,
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
        "U": 0,
    }

    # Test 2: CDS gene with unexpected amino acids in input dict
    test_gene_2 = CDSGene({"type": "cds", "id": "2", "gene": "test_gene_2", "product": "test_gene_product_2", "nt": "ATGCCGATTTAG", "aa": "MPIB"})
    test_gene_2.calculate_aa_counts()
    assert test_gene_2.aa_dict == {
        "C": 0,
        "Y": 0,
        "W": 0,
        "S": 0,
        "R": 0,
        "N": 0,
        "F": 0,
        "V": 0,
        "M": 1,
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
        "U": 0,
    }
    assert test_gene_2.unexpected_aas == ["B"]

def test_CDSGene_calculate_GCn():
    """Test code to calculate GCn content for a given gene."""
    test_gene = CDSGene({"type": "cds", "id": "1", "gene": "test_gene_1", "product": "test_gene_product_1", "nt": "CTGAATCGAACT", "aa": "LNRT"})
    test_gene.calculate_GCn()
    assert test_gene.GC1 == 0.5
    assert test_gene.GC2 == 0.5
    assert test_gene.GC3 == 0.25

def test_tRNAGene_init():
    """Test code to create tRNAGene object and initialize class attributes."""
    # Test 1: tRNA gene with all expected inputs
    test_gene = tRNAGene({"type": "tRNA", "id": "1", "gene": "test_gene_1", "product": "tRNA-Thr(cgt)", "score": 85.1, "nt": "GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA", "amino_acid": "Thr", "anti_codon": "cgt"})
    assert test_gene.input_error is False
    assert test_gene.type == "tRNA"
    assert test_gene.id == "1"
    assert test_gene.no_score is False
    assert test_gene.pseudogene is False
    assert test_gene.no_aa is False
    assert test_gene.no_anticodon is False
    assert not hasattr(test_gene, "unexpected_aa")
    assert not hasattr(test_gene, "unexpected_anticodon")
    assert test_gene.score == 85.1
    assert test_gene.amino_acid == "Thr"
    assert test_gene.anti_codon == "CGT"

    # Test 2: tRNA gene with missing score
    test_gene_2 = tRNAGene({"type": "tRNA", "id": "2", "gene": "test_gene_2", "product": "tRNA-Thr(cgt)", "nt": "GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA", "amino_acid": "Thr", "anti_codon": "cgt"})
    assert test_gene_2.input_error is False
    assert test_gene_2.type == "tRNA"
    assert test_gene_2.id == "2"
    assert test_gene_2.no_score is True

    # Test 3: tRNA gene that is a pseudogene
    test_gene_3 = tRNAGene({"type": "tRNA", "id": "3", "gene": "test_gene_3", "product": "tRNA-Xxx", "score": 21, "pseudogene": True, "nt": "GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA"})
    assert test_gene_3.input_error is False
    assert test_gene_3.type == "tRNA"
    assert test_gene_3.id == "3"
    assert test_gene_3.no_score is False
    assert test_gene_3.pseudogene is True

    # Test 4: tRNA gene with unexpected amino acid and anticodon
    test_gene_4 = tRNAGene({"type": "tRNA", "id": "4", "gene": "test_gene_4", "product": "tRNA-Xxx", "score": 21, "nt": "GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA", "amino_acid": "Xxx", "anti_codon": "nnn"})
    assert test_gene_4.input_error is False
    assert test_gene_4.type == "tRNA"
    assert test_gene_4.id == "4"
    assert test_gene_4.no_score is False
    assert test_gene_4.unexpected_aa == "Xxx"
    assert test_gene_4.unexpected_anti_codon == "nnn"
    assert not hasattr(test_gene_4, "amino_acid")
    assert not hasattr(test_gene_4, "anti_codon")


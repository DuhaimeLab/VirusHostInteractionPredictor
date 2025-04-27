"""Pytest for GeneSet methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import GeneSet


def test_GeneSet_init():
    """Test code to create GeneSet object and initialize class attributes."""
    # test 1 - test GeneSet object creation generates correct id and expected number of genes in genes attribute
    test_GeneSet = GeneSet("tests/datatests/test_annotated_genes.json")
    assert test_GeneSet.id == "test_annotated_genes"
    assert len(test_GeneSet.cds_genes) == 1
    assert len(test_GeneSet.tRNA_genes) == 1

    # test 2 - test Exception is raised if gene file is empty
    with pytest.raises(Exception):
        GeneSet("tests/datatests/test_empty_file.json")

    # test 3 - test Exception is raised if non-json file is provided
    with pytest.raises(Exception):
        GeneSet("tests/datatests/test_annotated_genes.ffn")

    # test 4 - test GeneSet object creation generates expected gene attributes for a good CDS gene
    assert test_GeneSet.cds_genes[0].input_error is False
    assert test_GeneSet.cds_genes[0].nt_input_error is False
    assert test_GeneSet.cds_genes[0].cds_len_error is False
    assert test_GeneSet.cds_genes[0].aa_input_error is False
    assert test_GeneSet.cds_genes[0].type == "cds"
    assert test_GeneSet.cds_genes[0].id == "good_CDS"
    assert test_GeneSet.cds_genes[0].gene == "good_CDS"
    assert test_GeneSet.cds_genes[0].product == "Chromosomal replication initiator protein DnaA"
    assert test_GeneSet.cds_genes[0].nt == "GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA"
    assert test_GeneSet.cds_genes[0].aa == "VSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS"
    assert test_GeneSet.cds_genes[0].codon_length == 3

    # test 5 - test GeneSet object creation generates expected gene attributes for a good tRNA gene
    assert test_GeneSet.tRNA_genes[0].input_error is False
    assert test_GeneSet.tRNA_genes[0].type == "tRNA"
    assert test_GeneSet.tRNA_genes[0].id == "good_tRNA"
    assert test_GeneSet.tRNA_genes[0].gene == "good_tRNA"
    assert test_GeneSet.tRNA_genes[0].product == "tRNA-Met(tca)"

    # test 6 - test genes are skipped correctly
    assert len(test_GeneSet.cds_general_input_errors) == 1
    assert len(test_GeneSet.cds_nt_input_errors) == 1
    assert len(test_GeneSet.cds_len_errors) == 1
    assert len(test_GeneSet.cds_aa_input_errors) == 1
    assert len(test_GeneSet.tRNA_general_input_errors) == 1
    assert len(test_GeneSet.tRNA_no_score) == 1
    assert len(test_GeneSet.tRNA_pseudogenes) == 1
    assert len(test_GeneSet.tRNA_no_aa) == 1
    assert len(test_GeneSet.tRNA_no_anticodon) == 1
    assert len(test_GeneSet.tRNA_unexpected_aa) == 1
    assert len(test_GeneSet.tRNA_unexpected_anti_codon) == 1
    assert test_GeneSet.skipped_cds_genes == 4/5
    assert test_GeneSet.skipped_tRNA_genes == 7/8


def test_GeneSet_codon_counts():
    """Test code to calculate codon counts across all CDS genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes.json")

    # test 1 - test high threshold_imprecise (100)
    test_GeneSet.codon_counts(threshold_imprecise = 100)
    assert len(test_GeneSet.imprecise_codons) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 0 # 100% of the codons would have to be unexpected for the gene to be skipped
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
        "sTGA": 0,
        "TGG": 0,
    }
    expected_aa_encoding_keys = set(test_GeneSet.codon_dict.keys()) - {"TGA", "TAA", "TAG"}
    assert set(test_GeneSet.codon_dict_aa_encoding.keys()) == expected_aa_encoding_keys

    # test 2 - test low threshold_imprecise (default, 0)
    test_GeneSet.codon_counts(threshold_imprecise = 0)
    assert len(test_GeneSet.imprecise_codons) == 1
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
        "sTGA": 0,
        "TGG": 0,
    }
    expected_aa_encoding_keys = set(test_GeneSet.codon_dict.keys()) - {"TGA", "TAA", "TAG"}
    assert set(test_GeneSet.codon_dict_aa_encoding.keys()) == expected_aa_encoding_keys


def test_GeneSet_codon_frequency():
    """Test code to calculate codon frequency across all CDS genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes.json")

    test_GeneSet.codon_frequency()  # using default threshold for tolerable imprecise codons
    assert len(test_GeneSet.imprecise_codons) == 1
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
        "sTGA": 0.0,
        "TGG": 0.0,
    }
    expected_aa_encoding_keys = set(test_GeneSet.codon_frq.keys()) - {"TGA", "TAA", "TAG"}
    assert set(test_GeneSet.codon_frq_aa_encoding.keys()) == expected_aa_encoding_keys


def test_GeneSet_amino_acid_counts():
    """Test code to calculate amino acid counts across all CDS genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes.json")

    # test 1 - test high threshold_unexpected (100)
    test_GeneSet.amino_acid_counts(threshold_unexpected = 100)
    assert len(test_GeneSet.unexpected_aas) == 1
    assert len(test_GeneSet.skipped_unexpected_peptides) == 0 # 100% of the amino acids would have to be unexpected for the gene to be skipped
    assert test_GeneSet.aa_dict == {
        "I": 0,
        "M": 1,
        "T": 0,
        "N": 0,
        "K": 0,
        "S": 3,
        "R": 0,
        "L": 1,
        "P": 0,
        "H": 0,
        "Q": 0,
        "V": 0,
        "A": 0,
        "D": 0,
        "E": 2,
        "G": 0,
        "F": 0,
        "Y": 0,
        "C": 0,
        "W": 0,
        "U": 0
    }

    # test 2 - test low threshold_unexpected (default, 0)
    test_GeneSet.amino_acid_counts(threshold_unexpected = 0)
    assert len(test_GeneSet.unexpected_aas) == 1
    assert len(test_GeneSet.skipped_unexpected_peptides) == 1
    assert test_GeneSet.aa_dict =={
        "I": 0,
        "M": 0,
        "T": 0,
        "N": 0,
        "K": 0,
        "S": 1,
        "R": 0,
        "L": 1,
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
        "U": 0
    }

def test_GeneSet_amino_acid_frequency():
    """Test code to calculate amino acid frequency across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes.json")

    test_GeneSet.amino_acid_frequency()  # using default threshold for tolerable unexpected aas
    assert len(test_GeneSet.unexpected_aas) == 1
    assert len(test_GeneSet.skipped_unexpected_peptides) == 1
    assert test_GeneSet.aa_frq == {
        "I": 0.0,
        "M": 0.0,
        "T": 0.0,
        "N": 0.0,
        "K": 0.0,
        "S": 1 / 3,
        "R": 0.0,
        "L": 1 / 3,
        "P": 0.0,
        "H": 0.0,
        "Q": 0.0,
        "V": 0.0,
        "A": 0.0,
        "D": 0.0,
        "E": 1 / 3,
        "G": 0.0,
        "F": 0.0,
        "Y": 0.0,
        "C": 0.0,
        "W": 0.0,
        "U": 0.0
    }


def test_GeneSet_RSCU():
    """Test code to calculate relative synonymous codon usage (RSCU) across all genes in a GeneSet object."""
    test_GeneSet = GeneSet("tests/datatests/test_short_genes.json")

    test_GeneSet.RSCU()  # using default threshold for tolerable imprecise codons
    assert len(test_GeneSet.imprecise_codons) == 1
    assert len(test_GeneSet.skipped_imprecise_genes) == 1
    assert test_GeneSet.RSCU_dict == {
        "ATA": 0.0,
        "ATC": 0.0,
        "ATT": 0.0,
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
    }

def test_GeneSet_tRNA_counts():
    """Test code to calculate tRNA gene copy counts across a GeneSet."""
    test_GeneSet = GeneSet("tests/datatests/test_tRNA_genes.json")
    test_GeneSet.tRNA_counts()
    assert test_GeneSet.tRNA_dict_aa == {
        "Ala": 0,
        "Arg": 2,
        "Asn": 0,
        "Asp": 0,
        "Cys": 0,
        "Gln": 0,
        "Glu": 0,
        "Gly": 0,
        "His": 0,
        "Ile": 0,
        "Leu": 0,
        "Lys": 0,
        "Met": 0,
        "Phe": 0,
        "Pro": 0,
        "Ser": 0,
        "Thr": 0,
        "Trp": 0,
        "Tyr": 0,
        "Val": 0,
        "SeC": 1,
        "fMet": 1,
        "Ile2": 1,
    }

    assert test_GeneSet.tRNA_dict_aa_1_letter == {
        "A": 0,
        "R": 2,
        "N": 0,
        "D": 0,
        "C": 0,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 0,
        "I": 1,
        "L": 0,
        "K": 0,
        "M": 1,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
        "U": 1,
    }

    assert test_GeneSet.tRNA_dict_tcc == {
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
    "AAA": 0,
    "AAG": 0,
    "AGC": 0,
    "AGT": 0,
    "AGA": 0,
    "AGG": 1,
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
    "CGT": 1,
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
    "TGC": 0,
    "TGT": 0,
    "sTGA": 1,
    "TGG": 0,
}


def test_GeneSet_tRNA_frequency():
    """Test code to calculate tRNA gene copy frequency out of total tRNA genes in a GeneSet."""
    test_GeneSet = GeneSet("tests/datatests/test_tRNA_genes.json")
    test_GeneSet.tRNA_frequency()
    assert test_GeneSet.tRNA_frq_aa == {
        "Ala": 0.0,
        "Arg": 2 / 5,
        "Asn": 0.0,
        "Asp": 0.0,
        "Cys": 0.0,
        "Gln": 0.0,
        "Glu": 0.0,
        "Gly": 0.0,
        "His": 0.0,
        "Ile": 0.0,
        "Leu": 0.0,
        "Lys": 0.0,
        "Met": 0.0,
        "Phe": 0.0,
        "Pro": 0.0,
        "Ser": 0.0,
        "Thr": 0.0,
        "Trp": 0.0,
        "Tyr": 0.0,
        "Val": 0.0,
        "SeC": 1 / 5,
        "fMet": 1 / 5,
        "Ile2": 1 / 5,
    }

    assert test_GeneSet.tRNA_frq_aa_1_letter == {
        "A": 0.0,
        "R": 2 / 5,
        "N": 0.0,
        "D": 0.0,
        "C": 0.0,
        "Q": 0.0,
        "E": 0.0,
        "G": 0.0,
        "H": 0.0,
        "I": 1 / 5,
        "L": 0.0,
        "K": 0.0,
        "M": 1 / 5,
        "F": 0.0,
        "P": 0.0,
        "S": 0.0,
        "T": 0.0,
        "W": 0.0,
        "Y": 0.0,
        "V": 0.0,
        "U": 1 / 5,
    }

    assert test_GeneSet.tRNA_frq_tcc == {
        "ATA": 0.0,
        "ATC": 0.0,
        "ATT": 0.0,
        "ATG": 2 / 5,
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
        "AGG": 1 / 5,
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
        "CGT": 1 / 5,
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
        "GAA": 0.0,
        "GAG": 0.0,
        "GGA": 0.0,
        "GGC": 0.0,
        "GGG": 0.0,
        "GGT": 0.0,
        "TCA": 0.0,
        "TCC": 0.0,
        "TCG": 0.0,
        "TCT": 0.0,
        "TTC": 0.0,
        "TTT": 0.0,
        "TTA": 0.0,
        "TTG": 0.0,
        "TAC": 0.0,
        "TAT": 0.0,
        "TGC": 0.0,
        "TGT": 0.0,
        "sTGA": 1 / 5,
        "TGG": 0.0,
    }

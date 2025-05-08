"""Test suite for the ComputeFeatures module."""

import math

import pandas as pd

from vhip.mlmodel.compute_ml_features import ComputeFeatures, Pairs
from vhip.mlmodel.gene_features import CodonBiasComparison, tRNAMetrics
from vhip.mlmodel.genomes_features import KmerProfile

test_virus_genome_dir = "tests/datatests/sequences/virus_genomes/"
test_host_genome_dir = "tests/datatests/sequences/host_genomes/"
test_virus_gene_dir = "tests/datatests/sequences/virus_genes/"
test_host_gene_dir = "tests/datatests/sequences/host_genes/"


def test_ComputeFeatures_list_genome_files():
    """Test that listed genome files are correct."""
    all_genome_filenames = [
        "GCA_003931015.1_ASM393101v1_genomic.fasta",
        "GCA_003927235.1_ASM392723v1_genomic.fasta",
        "GCA_003931415.1_ASM393141v1_genomic.fasta",
        "GCA_003344205.1_ASM334420v1_genomic.fasta",
        "GCA_005146815.1_ASM514681v1_genomic.fna.fasta",
        "GCA_001974575.1_ASM197457v1_genomic.fna.fasta",
        "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
    ]
    all_genomes = [
        "GCA_003931015.1_ASM393101v1_genomic",
        "GCA_003927235.1_ASM392723v1_genomic",
        "GCA_003931415.1_ASM393141v1_genomic",
        "GCA_003344205.1_ASM334420v1_genomic",
        "GCA_005146815.1_ASM514681v1_genomic.fna",
        "GCA_001974575.1_ASM197457v1_genomic.fna",
        "GCA_002875995.1_ASM287599v1_genomic.fna",
    ]

    all_genome_filenames.sort()
    all_genomes.sort()
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.all_genome_files.sort()
    test.all_genomes.sort()
    assert len(test.all_genome_files) == 7
    assert len(test.all_genomes) == 7
    assert len(test.virus_genome_filenames) == 4
    assert len(test.virus_genomes) == 4
    assert len(test.host_genome_filenames) == 3
    assert len(test.host_genomes) == 3
    assert test.all_genome_files == all_genome_filenames
    assert test.all_genomes == all_genomes


def test_ComputeFeatures_list_gene_files():
    """Test that listed gene files are correct."""
    all_gene_filenames = [
        "GCA_003931015.1_ASM393101v1_genomic.json",
        "GCA_003927235.1_ASM392723v1_genomic.json",
        "GCA_003931415.1_ASM393141v1_genomic.json",
        "GCA_003344205.1_ASM334420v1_genomic.json",
        "GCA_005146815.1_ASM514681v1_genomic.fna.json",
        "GCA_001974575.1_ASM197457v1_genomic.fna.json",
        "GCA_002875995.1_ASM287599v1_genomic.fna.json",
        "test_short_genes_file.json",
    ]
    all_genes = [
        "GCA_003931015.1_ASM393101v1_genomic",
        "GCA_003927235.1_ASM392723v1_genomic",
        "GCA_003931415.1_ASM393141v1_genomic",
        "GCA_003344205.1_ASM334420v1_genomic",
        "GCA_005146815.1_ASM514681v1_genomic.fna",
        "GCA_001974575.1_ASM197457v1_genomic.fna",
        "GCA_002875995.1_ASM287599v1_genomic.fna",
        "test_short_genes_file",
    ]

    all_gene_filenames.sort()
    all_genes.sort()
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.all_gene_files.sort()
    test.all_genes.sort()
    assert len(test.all_gene_files) == 8
    assert len(test.all_genes) == 8
    assert len(test.virus_gene_filenames) == 4
    assert len(test.virus_genes) == 4
    assert len(test.host_gene_filenames) == 4
    assert len(test.host_genes) == 4
    assert test.all_gene_files == all_gene_filenames
    assert test.all_genes == all_genes


def test_ComputeFeatures_pairs():
    """Test all pairs are generated correctly."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.determine_pairs()
    assert len(test.pairs) == 12
    assert all(isinstance(pair, Pairs) for pair in test.pairs)


def test_ComputeFeatures_get_headers():
    """Test that headers are retrieved correctly."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.get_headers()

    assert test.headers["MG592522.1"] == "GCA_003344205.1_ASM334420v1_genomic.fasta"
    assert (
        test.headers["MKKP01000001.1"]
        == "GCA_001974575.1_ASM197457v1_genomic.fna.fasta"
    )
    assert (
        test.headers["MKKP01000002.1"]
        == "GCA_001974575.1_ASM197457v1_genomic.fna.fasta"
    )
    assert (
        test.headers["MKKP01000003.1"]
        == "GCA_001974575.1_ASM197457v1_genomic.fna.fasta"
    )


def test_ComputeFeatures_add_blastn_files():
    """Test that blastn files are added correctly."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    blastn_filename = "tests/datatests/blastn.tsv"
    spacer_filename = "tests/datatests/spacer.tsv"
    test.add_blastn_files(blastn_filename, spacer_filename)
    assert test.blastn_path == blastn_filename
    assert test.spacer_path == spacer_filename


def test_ComputeFeatures_process_blastn():
    """Test that blastn file is processed correctly."""
    expected_results = {
        "GCA_003344205.1_ASM334420v1_genomic.fasta": [
            "GCA_003931015.1_ASM393101v1_genomic.fasta"
        ],
        "GCA_003931015.1_ASM393101v1_genomic.fasta": [
            "GCA_005146815.1_ASM514681v1_genomic.fna.fasta",
            "GCA_003931015.1_ASM393101v1_genomic.fasta",
        ],
    }
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.get_headers()
    test.add_blastn_files("tests/datatests/blastn_phagevhost.tsv", "")
    test.process_blastn()
    assert test.blastn == expected_results


def test_ComputeFeatures_process_spacers():
    """Test that spacers are processed correctly."""
    expected_results = {
        "GCA_003344205.1_ASM334420v1_genomic.fasta": [
            "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
            "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
            "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
        ],
        "GCA_003931015.1_ASM393101v1_genomic.fasta": [
            "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
            "GCA_002875995.1_ASM287599v1_genomic.fna.fasta",
        ],
    }
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.get_headers()
    test.add_blastn_files("", "tests/datatests/blastn_phagevspacer.tsv")
    test.process_spacers()
    assert test.spacers == expected_results


def test_ComputeFeatures_generate_kmer_profiles():
    """Test the kmer profiles are properly generated."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.determine_pairs()
    test.generate_kmer_profiles()
    assert isinstance(test.k6profiles, dict)
    assert isinstance(
        test.k6profiles["GCA_003344205.1_ASM334420v1_genomic"], KmerProfile
    )


def test_ComputeFeatures_generate_codon_profiles():
    """Test the codon counts and frequencies are properly generated with expected values."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.generate_codon_profiles()
    assert isinstance(test.codon_counts, dict)
    assert isinstance(test.codon_frqs, dict)
    assert isinstance(
        test.codon_counts["GCA_003344205.1_ASM334420v1_genomic"], dict
    )
    assert isinstance(test.codon_frqs["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert test.codon_counts["test_short_genes_file"] == {
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
    assert test.codon_frqs["test_short_genes_file"] == {
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
    expected_aa_encoding_keys = set(test.codon_counts["test_short_genes_file"].keys()) - {"TGA", "TAA", "TAG"}
    assert set(test.codon_counts_aa_encoding["test_short_genes_file"].keys()) == expected_aa_encoding_keys
    assert set(test.codon_frqs_aa_encoding["test_short_genes_file"].keys()) == expected_aa_encoding_keys

def test_ComputeFeatures_generate_aa_profiles():
    """Test the amino acid counts and frequencies are properly generated with expected values."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.generate_aa_profiles()
    assert isinstance(test.aa_counts, dict)
    assert isinstance(test.aa_frqs, dict)
    assert isinstance(
        test.aa_counts["GCA_003344205.1_ASM334420v1_genomic"], dict
    )
    assert isinstance(test.aa_frqs["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert test.aa_counts["test_short_genes_file"] == {
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
    assert test.aa_frqs["test_short_genes_file"] == {
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

def test_ComputeFeatures_generate_RSCU():
    """Test the relative synonymous codon usage (RSCU) profiles are properly generated with expected values."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.generate_RSCU()
    assert isinstance(test.RSCU, dict)
    assert isinstance(test.RSCU["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert test.RSCU["test_short_genes_file"] == {
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

def test_ComputeFeatures_determine_custom_pairs():
    """Test that custom pairs are generated correctly."""
    # test 1 - one custom pair
    test_one_pair = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
        pairs_of_interest="tests/datatests/custom_pairs_1.csv",
    )
    test_one_pair.list_genome_files()
    if test_one_pair.pairs_of_interest:
        test_one_pair.determine_custom_pairs(test_one_pair.pairs_of_interest)
    assert len(test_one_pair.pairs) == 1
    assert all(isinstance(pair, Pairs) for pair in test_one_pair.pairs)
    assert test_one_pair.pairs[0].virus == "GCA_003344205.1_ASM334420v1_genomic"
    assert (
        test_one_pair.pairs[0].host == "GCA_001974575.1_ASM197457v1_genomic.fna"
    )

    # test 2 - multiple custom pairs
    test_two_pairs = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
        pairs_of_interest="tests/datatests/custom_pairs_2.csv",
    )
    test_two_pairs.list_genome_files()
    if test_two_pairs.pairs_of_interest:
        test_two_pairs.determine_custom_pairs(test_two_pairs.pairs_of_interest)
    assert len(test_two_pairs.pairs) == 2
    assert all(isinstance(pair, Pairs) for pair in test_two_pairs.pairs)
    assert test_two_pairs.pairs[1].virus == "GCA_003344205.1_ASM334420v1_genomic"
    assert (
        test_two_pairs.pairs[1].host == "GCA_002875995.1_ASM287599v1_genomic.fna"
    )

    # test 3 - unavailable pair (missing virus or host)
    test_bad_pair = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
        pairs_of_interest="tests/datatests/custom_pairs_3.csv",
    )
    test_bad_pair.list_genome_files()
    if test_bad_pair.pairs_of_interest:
        test_bad_pair.determine_custom_pairs(test_bad_pair.pairs_of_interest)
    assert len(test_bad_pair.pairs) == 0

def test_ComputeFeatures_generate_tRNA_profiles():
    """Test the tRNA profiles are properly generated."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.generate_tRNA_profiles()
    assert isinstance(test.tRNA_counts_aa, dict)
    assert isinstance(test.tRNA_counts_aa_1_letter, dict)
    assert isinstance(test.tRNA_frqs_aa, dict)
    assert isinstance(test.tRNA_counts_tcc, dict)
    assert isinstance(test.tRNA_frqs_tcc, dict)
    assert isinstance(
        test.tRNA_counts_aa["GCA_003344205.1_ASM334420v1_genomic"], dict
    )
    assert isinstance(test.tRNA_counts_aa_1_letter["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert isinstance(test.tRNA_frqs_aa["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert isinstance(test.tRNA_frqs_aa_1_letter["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert isinstance(test.tRNA_counts_tcc["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert isinstance(test.tRNA_frqs_tcc["GCA_003344205.1_ASM334420v1_genomic"], dict)
    assert test.tRNA_counts_aa["test_short_genes_file"] == {
        "Ala": 0,
        "Arg": 0,
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
        "Ser": 1,
        "Thr": 0,
        "Trp": 0,
        "Tyr": 0,
        "Val": 0,
        "SeC": 0,
        "fMet": 0,
        "Ile2": 0,
    }

    assert test.tRNA_counts_aa_1_letter["test_short_genes_file"] == {
        "A": 0,
        "R": 0,
        "N": 0,
        "D": 0,
        "C": 0,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 0,
        "I": 0,
        "L": 0,
        "K": 0,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 1,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
        "U": 0,
    }

    assert test.tRNA_frqs_aa["test_short_genes_file"] == {
        "Ala": 0.0,
        "Arg": 0.0,
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
        "Ser": 1 / 1,
        "Thr": 0.0,
        "Trp": 0.0,
        "Tyr": 0.0,
        "Val": 0.0,
        "SeC": 0.0,
        "fMet": 0.0,
        "Ile2": 0.0,
    }

    assert test.tRNA_frqs_aa_1_letter["test_short_genes_file"] == {
        "A": 0.0,
        "R": 0.0,
        "N": 0.0,
        "D": 0.0,
        "C": 0.0,
        "Q": 0.0,
        "E": 0.0,
        "G": 0.0,
        "H": 0.0,
        "I": 0.0,
        "L": 0.0,
        "K": 0.0,
        "M": 0.0,
        "F": 0.0,
        "P": 0.0,
        "S": 1 / 1,
        "T": 0.0,
        "W": 0.0,
        "Y": 0.0,
        "V": 0.0,
        "U": 0.0,
    }

    assert test.tRNA_counts_tcc["test_short_genes_file"] == {
        "ATA": 0,
        "ATC": 0,
        "ATT": 0,
        "ATG": 0,
        "ACA": 0,
        "ACC": 0,
        "ACG": 0,
        "ACT": 0,
        "AAC": 0,
        "AAT": 0,
        "AAA": 0,
        "AAG": 0,
        "AGC": 1,
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
        "sTGA": 0,
        "TGG": 0,
    }

    assert test.tRNA_counts_tcc["test_short_genes_file"] == {
        "ATA": 0.0,
        "ATC": 0.0,
        "ATT": 0.0,
        "ATG": 0.0,
        "ACA": 0.0,
        "ACC": 0.0,
        "ACG": 0.0,
        "ACT": 0.0,
        "AAC": 0.0,
        "AAT": 0.0,
        "AAA": 0.0,
        "AAG": 0.0,
        "AGC": 1 / 1,
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
        "sTGA": 0.0,
        "TGG": 0.0,
    }


def test_ComputeFeatures_compute_features():
    """Test all pair properties are populated correctly from running compute_feature()."""
    test_CF = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
        pairs_of_interest="tests/datatests/custom_pairs_1.csv",
    )
    test_CF.add_blastn_files(
        "tests/datatests/blastn_phagevhost.tsv",
        "tests/datatests/blastn_phagevspacer.tsv",
    )
    test_CF.do_setup()  # determining virus-host pairs (in this case, initialize one custom pair), gets fasta headers, read and process blastn_output, compute GC content and k-mer profiles, and for each organism: generate dictionaries of codon, amino acid, and synonymous codon usage frequencies.
    test_CF.compute_features(
        test_CF.pairs[0]
    )  # computes comparisons (e.g. distances) between profiles generated in do_setup for this single virus and single host pair

    # Genome-level features
    assert math.isclose(
        test_CF.pairs[0].GCdifference, -1.7559344991201087, rel_tol=1e-6
    )
    assert math.isclose(
        test_CF.pairs[0].GCdifference, -1.7559344991201087, rel_tol=1e-6
    )
    assert math.isclose(test_CF.pairs[0].k3dist, 0.07392896302485052, rel_tol=1e-6)
    assert math.isclose(test_CF.pairs[0].k6dist, 0.18284169630469121, rel_tol=1e-6)
    assert test_CF.pairs[0].homology_hit is False

    # Gene-level features
    # Codon bias comparison
    assert isinstance(test_CF.pairs[0].codons_comparison, CodonBiasComparison)
    assert math.isclose(test_CF.pairs[0].codons_comparison.slope, 0.9718185281167805, rel_tol=1e-6)  #mystery pytest failure
    assert math.isclose(
        test_CF.pairs[0].codons_comparison.R2, 0.768361516064103, rel_tol=1e-6
    )
    assert math.isclose(
        test_CF.pairs[0].codons_comparison.cos_similarity,
        0.9632986150103935,
        rel_tol=1e-6,
    )

    # Amino acid bias comparison
    assert isinstance(test_CF.pairs[0].aa_comparison, CodonBiasComparison)
    assert math.isclose(test_CF.pairs[0].aa_comparison.slope, 0.9317915647675197, rel_tol=1e-6) #mystery pytest failure
    assert math.isclose(
        test_CF.pairs[0].aa_comparison.R2, 0.8959324233050872, rel_tol=1e-6
    )
    assert math.isclose(
        test_CF.pairs[0].aa_comparison.cos_similarity, 0.9888428425307265, rel_tol=1e-6
    )

    # RSCU comparison
    assert isinstance(test_CF.pairs[0].RSCU_comparison, CodonBiasComparison)
    assert math.isclose(test_CF.pairs[0].RSCU_comparison.slope, 0.5981061832691094, rel_tol=1e-6) #mystery pytest failure
    assert math.isclose(
        test_CF.pairs[0].RSCU_comparison.R2, 0.40950550601057134, rel_tol=1e-6
    )
    assert math.isclose(
        test_CF.pairs[0].RSCU_comparison.cos_similarity,
        0.9465383699397703,
        rel_tol=1e-6,
    )

    # TAAI
    assert isinstance(test_CF.pairs[0].tRNAMetrics, tRNAMetrics)
    assert math.isclose(test_CF.pairs[0].tRNAMetrics.virusTAAI_hosttRNA, 0.17806764388058136, rel_tol=1e-6)
    assert math.isclose(test_CF.pairs[0].tRNAMetrics.virusTAAI_totaltRNA, 0.17806764388058136, rel_tol=1e-6)

    # TCAI
    assert math.isclose(test_CF.pairs[0].tRNAMetrics.virusTCAI_hosttRNA, -0.06738491986617158, rel_tol=1e-6)
    assert math.isclose(test_CF.pairs[0].tRNAMetrics.virusTCAI_totaltRNA, -0.06738491986617158, rel_tol=1e-6)


def test_ComputeFeatures_complete_pipeline():
    """Check the complete pipeline for ComputeFeatures is working as intended."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.add_blastn_files(
        "tests/datatests/blastn_phagevhost.tsv",
        "tests/datatests/blastn_phagevspacer.tsv",
    )
    test.do_setup()
    test.run_parallel()
    test.convert_to_dataframe()
    assert isinstance(test.features_df, pd.DataFrame)
    assert test.features_df.shape[0] == 12
    assert test.features_df.shape[1] == 17

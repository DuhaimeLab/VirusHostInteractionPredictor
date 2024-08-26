"""Test suite for the ComputeFeatures module."""

import pandas as pd

from vhip.mlmodel.compute_ml_features import ComputeFeatures, Pairs
from vhip.mlmodel.genomes_features import KmerProfile
from vhip.mlmodel.gene_features import GeneSet

test_virus_genome_dir = "tests/datatests/sequences/virus_genomes/"
test_host_genome_dir = "tests/datatests/sequences/host_genomes/"
test_virus_gene_dir = "tests/datatests/sequences/virus_genes"
test_host_gene_dir = "tests/datatests/sequences/host_genes"


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

    all_genome_filenames.sort()
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.all_genome_files.sort()
    assert len(test.all_genome_files) == 7
    assert len(test.virus_genome_filenames) == 4
    assert len(test.host_genome_filenames) == 3
    assert test.all_genome_files == all_genome_filenames


def test_ComputeFeatures_list_gene_files():
    """Test that listed gene files are correct."""
    all_gene_filenames = [
        "GCA_003931015.1_ASM393101v1_genomic.ffn",
        "GCA_003927235.1_ASM392723v1_genomic.ffn",
        "GCA_003931415.1_ASM393141v1_genomic.ffn",
        "GCA_003344205.1_ASM334420v1_genomic.ffn",
        "GCA_005146815.1_ASM514681v1_genomic.fna.ffn",
        "GCA_001974575.1_ASM197457v1_genomic.fna.ffn",
        "GCA_002875995.1_ASM287599v1_genomic.fna.ffn",
    ]

    all_gene_filenames.sort()
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.all_gene_files.sort()
    assert len(test.all_gene_files) == 7
    assert len(test.virus_gene_filenames) == 4
    assert len(test.host_gene_filenames) == 3
    assert test.all_gene_files == all_gene_filenames


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
        test.k6profiles["GCA_003344205.1_ASM334420v1_genomic.fasta"], KmerProfile
    )

def test_ComputeFeatures_generate_codon_frq():
    """Test the codon frequency profiles are properly generated."""
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_gene_files()
    test.generate_codon_frq()
    assert isinstance(test.codon_frqs, dict)
    assert isinstance(test.codon_counts, dict)
    assert isinstance(
        test.codon_frqs["GCA_003344205.1_ASM334420v1_genomic.ffn"], dict
    )
    assert isinstance(
        test.codon_counts["GCA_003344205.1_ASM334420v1_genomic.ffn"], dict
    )

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
    assert test.features_df.shape[1] == 4


if __name__ == "__main__":
    test = ComputeFeatures(
        test_virus_genome_dir,
        test_host_genome_dir,
        test_virus_gene_dir,
        test_host_gene_dir,
    )
    test.list_genome_files()
    test.determine_pairs()
    test.generate_kmer_profiles()
    print(test.k6profiles["GCA_003344205.1_ASM334420v1_genomic.fasta"])

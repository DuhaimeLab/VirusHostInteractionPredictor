"""Contains classes to compute gene-level features.

This module provides:
- Gene: calculate codon counts, amino acid counts, and imprecise codon counts for a single gene
- GeneSet: calculate codon counts, codon frequency, RSCU (relative synonymous codon usage), amino acid counts, and amino acid frequency for a set of genes in an annotated file
- CodonBiasComparison: compare codon bias measurements between two gene sets using linear regression (slope, R^2) and cosine similarity
"""

import os
from typing import List

from .read_sequence import read_gene_products, read_headers, read_sequence

# Set up Codon Table with each codon's encoded amino acid (1 letter abbreviation)
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}

# Separate CODON_TABLE dictionary into lists of codons and amino acids
CODON_LIST = CODON_TABLE.keys()
AA_LIST = set(CODON_TABLE.values())


# Define Gene class
class Gene:
    """Class representing a gene.

    Args:
        gene_seq (str): The nucleotide sequence of the gene.
        codon_length (int): Length of 1 codon (default is 3).
    """

    def __init__(self, gene_seq: str, codon_length: int = 3, gene_id: str = '', gene_product: str = '') -> None:
        """Initialize class variables."""
        if len(gene_seq) % codon_length == 0:
            self.seq = gene_seq
            self.codon_length = codon_length
            self.gene_id = gene_id
            self.gene_product = gene_product
        elif len(gene_seq) % codon_length != 0:
            raise Exception("Gene length is not a multiple of codon length.")

    def calculate_codon_counts(self) -> None:
        """Calculate counts of each unique codon in a gene.

        Populates the following class attributes:
            self.codon_dict (str: int): Each key of dictionary is a unique codon, and the values represent the number of times the associated codon (key) appears in the provided gene sequence.
            self.number_imprecise_codons (int): Number of codons that are not precise (i.e. are not found in expected CODON_LIST).
        """
        self.number_imprecise_codons = 0
        self.codon_dict = dict.fromkeys(CODON_LIST, 0)

        for i in range(0, len(self.seq), self.codon_length):
            codon = self.seq[i : i + self.codon_length]
            if codon in self.codon_dict.keys():
                self.codon_dict[codon] += 1
            else:
                self.number_imprecise_codons += 1

    def calculate_aa_counts(self) -> None:
        """Calculate counts of each unique amino acid encoded by a gene.

        Populates the following class attributes:
            self.aa_dict (str: int): Each key of dictionary is an unique amino acid, and values represent the number of times the associated amino acid (key) appears to be encoded by codons in the gene sequence.
        If not previously calculated, method will also populate:
            self.codon_dict (str: int): Each key of dictionary is a unique codon, and the values represent the number of times the associated codon (key) appears in the provided gene sequence.
            self.number_imprecise_codons (int): Number of codons that are not precise (i.e. are not found in expected CODON_LIST).
        """
        self.aa_dict = dict.fromkeys(AA_LIST, 0)

        if not hasattr(self, "codon_dict"):
            self.calculate_codon_counts()

        for codon in self.codon_dict:
            if self.codon_dict[codon] != 0:
                aa = CODON_TABLE[codon]
                self.aa_dict[aa] += self.codon_dict[codon]


# Define GeneSet class
class GeneSet:
    """Class representing a gene set, usually the genes predicted from a genome sequence.

    Args:
        gene_file (str): Path of annotated genes file containing gene set of interest.
    """

    def __init__(self, gene_file: str) -> None:
        """Initialize class variables and read in an annotated genes file, storing Gene objects and metadata in lists."""
        if not gene_file or not os.path.getsize(gene_file):
            raise Exception(
                "Genes file is not provided or empty. Please provide a valid gene file."
            )

        sequences = read_sequence(gene_file)
        self.genes: List[Gene] = [Gene(seq) for seq in sequences]
        self.gene_ids: List[str] = read_headers(gene_file)
        self.gene_products: List[str] = read_gene_products(gene_file)

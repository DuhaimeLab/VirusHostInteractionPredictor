"""Contains classes to compute gene-level features.

This module provides:
- Gene: calculate codon counts, amino acid counts, and imprecise codon counts for a single gene
- GeneSet: calculate codon counts, codon frequency, RSCU (relative synonymous codon usage), amino acid counts, and amino acid frequency for a set of genes in an annotated file
- CodonBiasComparison: compare codon bias measurements between two gene sets using linear regression (slope, R^2) and cosine similarity
"""

import os
import re
from typing import List, Union

import numpy as np
import scipy  # pyright: ignore[reportMissingTypeStubs]

from .read_sequence import read_annotated_genes, reverse_complement

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

# Separate CODON_TABLE dictionary into lists of codons, amino acids, stop codons, and non-degenerate codons (encoded amino acid is specific to one codon alone)
CODON_LIST = list(CODON_TABLE.keys())
AA_LIST = [aa for aa in CODON_TABLE.values() if aa != "_"]
stop_codons = [codon for codon, aa in CODON_TABLE.items() if aa == "_"]
non_degenerate_codons = [
    codon
    for codon, aa in CODON_TABLE.items()
    if list(CODON_TABLE.values()).count(aa) == 1
]

# Amino acid abbreviations conversions
AA_CONVERSIONS = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}


# Define Gene class
class Gene:
    """Class representing a gene.

    Args:
        json_dict (dict): Dictionary containing gene information. Required keys are 'type', 'id', 'gene', 'product'. Descriptions below:
            - type (str): Type of gene ('cds' and 'tRNA' values will inform annotation parsing).
            - id (str): Gene ID.
            - gene (str): Gene name.
            - product (str): Gene product name.

    Populates the following class attributes:
        self.input_error (bool): True if input dictionary does not contain all expected keys, wherein method will exit (expect False).
        self.type (str),
        self.id (str),
        self.gene (str),
        self.product (str),
    """
    def __init__(
        self,
        json_dict: dict[str, str]
    ) -> None:
        """Initialize class variables."""
        self.input_error: bool = False # to flag any missing generic keys from json input dict
        if not all(key in json_dict.keys() for key in ["type", "id", "gene", "product"]): # confirm all required keys present in input dictionary
            print("Input dictionary does not contain 'type', 'id', 'gene', 'product' keys. See documentation for Gene class initialization.")
            self.input_error = True
            return
        else: # populate basic class attributes if provided
            self.type: str = json_dict["type"]
            self.id: str = json_dict["id"]
            self.gene: str = json_dict["gene"]
            self.product: str = json_dict["product"]


# Define CDSGene subclass
class CDSGene(Gene):
    """Class representing a CDS gene (sublass of Gene).

    Args:
        json_dict (dict): Dictionary containing gene information. Required keys are 'type', 'id', 'gene', 'product', 'nt', and 'aa'. Descriptions below:
            - type (str): Type of gene (must be 'cds', otherwise error).
            - id (str): Gene ID.
            - gene (str): Gene name.
            - product (str): Gene product name.
            - nt (str): Nucleotide sequence of the gene. This string should be divisible by the codon length (default 3). Note that degenerate codons will be skipped.
            - aa (str): Amino acid sequence of the gene.

    Populates the following class attributes:
        self.input_error (bool): True if input dictionary does not contain expected general info keys ('type', 'id', 'gene', 'product'), wherein method will exit (expect False).
        self.nt_input_error (bool): True if the nucleotide sequence is not provided in the input dictionary, wherein method will exit (expect False).
        self.cds_len_error (bool): True if the length of the nucleotide sequence of a CDS is not divisible by the codon length, wherein method will exit (expect False).
        self.aa_input_error (bool): True if the amino acid sequence is not provided in the input dictionary, wherein method will exit (expect False).
        self.type (str),
        self.id (str),
        self.gene (str),
        self.product (str),
        self.nt (str),
        self.codon_length (int): Length of 1 codon (default is 3 - if user specifies an alternative, custom codon dictionary must be provided in source code to compute codon counts).
        self.n_codons (int): Number of codons in the nucleotide sequence.
        self.aa (str): Populated if a CDS gene,
    """
    def __init__(self, json_dict: dict[str, str]) -> None:
        """Initialize class variables."""
        self.codon_length: int = 3
        self.nt_input_error: bool = False # will flag if nucleotide sequence not provided in json input dict
        self.cds_len_error: bool = False # will flag if length of nucleotide sequence not divisible by codon length
        self.aa_input_error: bool = False # will flag if amino acid sequence not provided in json input dict
        super().__init__(json_dict)

        if self.input_error is True: # exit if input dictionary does not contain all required keys
            return
        elif json_dict["type"] != "cds": # exit if input gene type is not 'cds'
            raise Exception("Gene is not a CDS gene (expected dictionary element 'type': 'cds'). This class is for CDS genes only.")
        elif "nt" not in json_dict.keys(): # exit if nucleotide sequence not provided for CDS gene
            print(f"{json_dict["id"]}: Input dictionary does not contain 'nt' key for cds gene. See documentation for Gene class initialization.")
            self.nt_input_error: bool = True
            return
        elif len(json_dict["nt"]) % self.codon_length != 0: # exit if gene length not divisible by 3 for CDS gene
            print(f"{json_dict["id"]}: Length of nucleotide sequence is not divisible by codon length.")
            self.cds_len_error: bool = True
            return
        elif "aa" not in json_dict.keys(): # exit if amino acid sequence not provided for CDS gene
            print(f"{json_dict["id"]}: Input dictionary does not contain 'aa' key for cds gene. See documentation for Gene class initialization.")
            self.aa_input_error: bool = True
            return
        else: # populate CDS gene attributes if provided
            self.nt = json_dict["nt"]
            self.n_codons: int = len(self.nt) // self.codon_length
            self.aa: str = json_dict["aa"]

    def calculate_codon_counts(self) -> None:
        """Calculate counts of each unique codon in a cds gene.

        Populates the following class attributes:
            self.codon_dict (str: int): Each key of dictionary is a unique codon, and the values represent the number of times the associated codon (key) appears in the provided gene sequence.
            self.imprecise_codons (int): Percentage of codons that are not precise (i.e. are not found in expected CODON_LIST and may contain degeneracies).
        """
        self.codon_dict = dict.fromkeys(CODON_LIST, 0)
        self.imprecise_codons: List[str] = []

        if not (self.input_error or self.nt_input_error or self.cds_len_error or self.aa_input_error):
            for i in range(0, len(self.nt), self.codon_length):
                codon = self.nt[i : i + self.codon_length]
                if codon in self.codon_dict.keys():
                    self.codon_dict[codon] += 1
                else:
                    self.imprecise_codons.append(codon)

    def calculate_aa_counts(self) -> None:
        """Calculate counts of each unique amino acid encoded by a cds gene.

        Populates the following class attributes:
            self.aa_dict (str: int): Each key of dictionary is an unique amino acid, and values represent the number of times the associated amino acid (key) appears to be encoded by codons in the gene sequence.
        """
        self.aa_dict = dict.fromkeys(AA_LIST, 0)
        self.unexpected_aas: List[str] = []

        if not (self.input_error or self.nt_input_error or self.cds_len_error or self.aa_input_error):
            for aa in self.aa:
                if aa in self.aa_dict.keys():
                    self.aa_dict[aa] += 1
                elif aa not in self.aa_dict.keys():
                    self.unexpected_aas.append(aa)

    def calculate_GCn(self) -> None:
        """Calculate GC content at position 1, 2, and 3 of a gene.

        Populates the following class attributes:
            self.GC1 (float): GC content of the gene at position 1.
            self.GC2 (float): GC content of the gene at position 2.
            self.GC3 (float): GC content of the gene at position 3.
        """
        gc1: int = 0
        gc2: int = 0
        gc3: int = 0

        for i in range(0, len(self.nt), self.codon_length):
            codon = self.nt[i : i + self.codon_length]
            for j in range(self.codon_length):
                if codon[j] == "G" or codon[j] == "C":
                    if j == 0:
                        gc1 += 1
                    elif j == 1:
                        gc2 += 1
                    elif j == 2:
                        gc3 += 1

        if self.n_codons > 0:
            self.GC1 = gc1 / self.n_codons
            self.GC2 = gc2 / self.n_codons
            self.GC3 = gc3 / self.n_codons

# Define tRNAGene subclass
class tRNAGene(Gene):
    """Class representing a tRNA gene (sublass of Gene).

    Args:
        json_dict (dict): Dictionary containing gene information. Required keys are 'type', 'id', 'gene', 'product', 'score', 'amino_acid', and 'anti_codon'. Descriptions below:
            - type (str): Type of gene (must be 'cds', otherwise error).
            - id (str): Gene ID.
            - gene (str): Gene name.
            - product (str): Gene product name.
            - score (float): Bakta output confidence score for the tRNA gene prediction.
            - amino_acid (str): Required unless pseudogene. Amino acid associated with the tRNA gene. Only those in with an associated element in AA_CONVERSIONS will be counted in downstream methods.
            - anti_codon (str): Required unlesss psueodgene. Anticodon sequence associated with the tRNA gene.

    Populates the following class attributes:
        self.input_error (bool): True if input dictionary does not contain expected general Gene info keys ('type', 'id', 'gene', 'product'), wherein method will exit (expect False).
        self.no_score (bool): True if the score is not provided in the input dictionary, wherein method will exit (expect False).
        self.pseudogene (bool): True if the tRNA gene is a pseudogene (output listed as tRNA-Xxx with no amino_acid and anti_codon keys).
        self.no_aa (bool): True if the amino acid is not provided for a non-pseudogene, wherein method will exit (expect False).
        self.no_anticodon (bool): True if the anticodon is not provided for a non-pseudogene, wherein method will exit (expect False).
        self.unexpected_aa (List['str']): True if provided amino_acid value is not in AA_CONVERSIONS, wherein method will exit.
        self.unexpected_anti_codon (List['str']): True if anticodon value does not only contain 'a','t','g','c' substrings, wherein method will exit.
        self.type (str),
        self.id (str),
        self.gene (str),
        self.product (str),
        self.score (float),
        self.amino_acid (str),
        self.anti_codon (str)
    """
    def __init__(self, json_dict: dict[str, str]) -> None:
        """Initialize class variables."""
        self.no_score: bool = False # will flag if score not provided in json input dict
        self.pseudogene: bool = False # will flag if tRNA gene is a pseudogene
        self.no_aa: bool = False # will flag if amino acid value not provided for non-pseudogene
        self.no_anticodon: bool = False # will flag if anticodon value not provided for non-pseudogene
        self.unexpected_aa: bool = False # will flag if amino acid value is not in AA_CONVERSIONS
        self.unexpected_anti_codon: bool = False # will flag if anticodon value does not only contain 'a','t','g','c' substrings

        # Initialize superclass Gene attributes
        super().__init__(json_dict)

        # Check for required keys in input dictionary
        if self.input_error is True: # exit if input dictionary does not contain all required general Gene keys
            return
        elif json_dict["type"] != "tRNA": # exit if input gene type is not 'cds'
            raise Exception("Gene is not a tRNA gene (expected dictionary element 'type': 'tRNA'). This class is for CDS genes only.")
        elif "score" not in json_dict.keys(): # exit if score not provided for tRNA gene
            self.no_score = True
            print(f"{json_dict["id"]}: Input dictionary does not contain 'score' key for tRNA gene.")
            return
        else:
            self.score = json_dict["score"] # populate score attribute if provided

        # Identify pseudogenes and populate tRNA gene attributes for non-pseudogenes
        if "pseudogene" in json_dict.keys(): # flag pseudogenes
            self.pseudogene = True
            print(f"{json_dict["id"]}: tRNA gene is a pseudogene.")
        elif "pseudogene" not in json_dict.keys(): # populate tRNA gene attributes appropriate for non-pseudogenes, flagging unexpected inputs
            # Amino acid attributes:
            if "amino_acid" not in json_dict.keys():
                self.no_aa = True
                print(f"{json_dict["id"]}: Input dictionary does not contain 'amino_acid' key for tRNA gene.")
                return
            elif json_dict["amino_acid"] not in AA_CONVERSIONS.keys():
                self.unexpected_aa = True
                print(f"{json_dict["id"]}: Unexpected amino acid provided for tRNA gene.")
                return
            else:
                self.amino_acid = json_dict["amino_acid"]
            # Anticodon attributes:
            if "anti_codon" not in json_dict.keys():
                self.no_anticodon = True
                print(f"{json_dict["id"]}: Input dictionary does not contain 'anti_codon' key for tRNA gene.")
                return
            elif not all(base in "atgc" for base in json_dict["anti_codon"].lower()):
                self.unexpected_anti_codon = True
                print(f"{json_dict["id"]}: Unexpected anticodon provided for tRNA gene.")
                return
            else:
                self.anti_codon = json_dict["anti_codon"].upper()

# Define GeneSet class
class GeneSet:
    """Class representing a gene set, usually the genes predicted from a genome sequence.

    Args:
        gene_file (str): Path of annotated genes file containing gene set of interest. .json output format from bakta.
    """

    def __init__(self, gene_file: str) -> None:
        """Initialize class variables and read in an annotated genes file, storing Gene objects and metadata in lists."""
        # Raise exceptions if gene_file input is not expected
        if not gene_file or not os.path.getsize(gene_file):
            raise Exception(
                "Genes file is not provided or empty. Please provide a valid gene file."
            )
        if not gene_file.endswith(".json"):
            raise Exception(
                "Gene file is not in .json format. Please provide a valid gene file."
            )
        # Initialize class attributes
        self.id = os.path.splitext(os.path.basename(gene_file))[0]
        self.cds_genes: List[CDSGene] = []
        self.tRNA_genes: List[Gene] = []
        self.skipped_cds_genes: float = 0
        self.skipped_tRNA_genes: float = 0

        # Initialize Quality control attributes
        self.cds_general_input_errors: List[Gene] = []
        self.cds_nt_input_errors: List[Gene] = []
        self.cds_len_errors: List[Gene] = []
        self.cds_aa_input_errors: List[Gene] = []
        self.tRNA_input_errors: List[Gene] = []

        # Store tRNA/CDS genes, filtering out genes with unexpected or missing inputs
        all_input_genes = read_annotated_genes(gene_file)
        self.all_input_cds_genes = [gene for gene in all_input_genes if gene["type"] == "cds"]
        self.all_input_tRNA_genes = [gene for gene in all_input_genes if gene["type"] == "tRNA"]
        for gene in self.all_input_cds_genes:
            current_gene = CDSGene(gene)
            if current_gene.input_error is True:
                self.cds_general_input_errors.append(current_gene)
            elif current_gene.nt_input_error is True:
                self.cds_nt_input_errors.append(current_gene)
            elif current_gene.cds_len_error is True:
                self.cds_len_errors.append(current_gene)
            elif current_gene.aa_input_error is True:
                self.cds_aa_input_errors.append(current_gene)
            else:
                self.cds_genes.append(CDSGene(gene))
        for gene in self.all_input_tRNA_genes:
            current_gene = Gene(gene)
            if current_gene.input_error is True:
                self.tRNA_input_errors.append(current_gene)
            else:
                self.tRNA_genes.append(Gene(gene))

        # Report skipped genes as a fraction of GeneSet (if no input genes, skip attributes will not reflect this...)
        if len(self.all_input_cds_genes) > 0:
            n_skipped_cds_genes = len(self.cds_general_input_errors) + len(self.cds_nt_input_errors)+ len(self.cds_len_errors) + len(self.cds_aa_input_errors)
            self.skipped_cds_genes = n_skipped_cds_genes / len(self.all_input_cds_genes)
            print(f"{self.skipped_cds_genes * 100}% of CDS genes in {self.id} skipped due to missing info and/or lack of divisibility by codon length.")
        else:
            print(f"No input CDS genes in {self.id}.")
        if len(self.all_input_tRNA_genes) > 0:
            self.skipped_tRNA_genes = len(self.tRNA_input_errors) / len(self.all_input_tRNA_genes)
            print(f"{self.skipped_tRNA_genes * 100}% of tRNA genes in {self.id} skipped due to missing info.")
        else:
            print(f"No input tRNA genes in {self.id}.")

    def codon_counts(
        self, threshold_imprecise: float = 0.0
    ) -> None:
        """Aggregate the counts for each unique codon and imprecise codons across all CDS genes in a GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single CDS gene included in the GeneSet (default 0.0 or 0%)
        Populates the following class attributes:
            self.codon_dict (str: int): Counts of each unique codon across all CDS genes in the GeneSet.
            self.imprecise_codons (List(str)): list of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[Gene]): List of CDSGenes in the GeneSet that have more than threshold_imprecise codons.
        """
        # Check if GeneSet has any CDS genes
        if len(self.cds_genes) == 0:
            print(f"No Valid CDS genes in {self.id}. Skipping codon counting.")
            return

        # Initialize attributes
        self.codon_dict: dict[str, int] = dict.fromkeys(CODON_LIST, 0)
        self.imprecise_codons: List[str] = []
        self.skipped_imprecise_genes: List[Gene] = []

        counter = 0
        for gene in self.cds_genes:
            counter += 1
            if gene.n_codons < 1:
                print(f"{gene.id} of {self.id} has no codons. Skipping.")
                continue

            gene.calculate_codon_counts() # calculate codon counts for the current gene
            self.imprecise_codons.extend(gene.imprecise_codons)  # add all elements of gene.imprecise_codons to GeneSet imprecise codons

            # if percentage of codons in current gene is over threshold, skip the gene
            if len(gene.imprecise_codons)/gene.n_codons <= threshold_imprecise:
                for key, val in gene.codon_dict.items():
                    self.codon_dict[key] += val
            else:
                self.skipped_imprecise_genes.append(gene)

    def codon_frequency(
        self, threshold_imprecise: float = 0.0
    ) -> None:
        """Calculate the frequency of each unique codon across all CDS genes in a GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single CDS gene included in the GeneSet (default 0.0 or 0%)
        Populates the following class attributes:
            self.codon_frq (str: float): Frequency of each unique codon across all CDS genes in the GeneSet.
        If not populated previously by running codon_counts():
            self.codon_dict (str: int): Counts of each unique codon across all CDS genes in the GeneSet.
            self.imprecise_codons (List(str)): list of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[Gene]): List of CDSGenes in the GeneSet that have more than threshold_imprecise codons.
        """
        # Check if GeneSet has any CDS genes
        if len(self.cds_genes) == 0:
            print(f"No Valid CDS genes in {self.id}. Skipping codon frequency calculation.")
            return

        print(f"Calculating codon frequencies in {self.id}.")
        self.codon_frq: dict[str, float] = {}

        if not hasattr(self, "codon_dict"):
            # If aggregate codon counts have not already been calculated, runs codon_counts()
            self.codon_counts(
                threshold_imprecise=threshold_imprecise
            )

        if hasattr(self, "codon_dict") and len(self.skipped_imprecise_genes)/len(self.cds_genes) < 1:
            total = sum(self.codon_dict.values())
            if total >= 1:
                self.codon_frq = {k: (v / total) for k, v in self.codon_dict.items()}
            elif total < 1:
                print(f"No valid codons in {self.id}. Cannot calculate frequencies.")
        else:
            print(f"No valid and precise CDS genes in {self.id} to calculate codon frequencies.")

    def amino_acid_counts(
        self, threshold_unexpected: float = 0.0
    ) -> None:
        """Aggregate the counts for each unique amion acid and unexpected amino acid across all CDS genes in a GeneSet.

        Args:
            threshold_unexpected (float): Percentage of unexpected (not in AA_LIST) amino acids tolerated in a single CDS gene included in the GeneSet (default 0.0 or 0%)
        Populates the following class attributes:
            self.aa_dict (str: int): Counts of each unique amino acid across all CDS genes in the GeneSet.
            self.unexpected_aas (list(str)): list of unexpected amino acids found in the GeneSet.
            self.skipped_unexpected_peptides (List[Gene]): List of CDSGenes in the GeneSet that have more than threshold_unexpected amino acids.
        """
        # Check if GeneSet has any CDS genes
        if len(self.cds_genes) == 0:
            print(f"No Valid CDS genes in {self.id}. Skipping amino acid counting.")
            return

        # Initialize attributes
        self.aa_dict: dict[str, int] = dict.fromkeys(AA_LIST, 0)
        self.unexpected_aas: List[str] = []
        self.skipped_unexpected_peptides: List[Gene] = []

        counter = 0
        for gene in self.cds_genes:
            counter += 1
            if len(gene.aa) < 1:
                print(f"{gene.id} of {self.id} has no reported amino acids. Skipping.")
                continue

            gene.calculate_aa_counts() # calculate amino acid counts for the current gene
            self.unexpected_aas.extend(gene.unexpected_aas)  # add all elements of gene.imprecise_codons to GeneSet imprecise codons

            # if percentage of amion acids in current gene is over threshold, skip the gene
            if len(gene.unexpected_aas)/len(gene.aa) <= threshold_unexpected:
                for key, val in gene.aa_dict.items():
                    self.aa_dict[key] += val
            else:
                self.skipped_unexpected_peptides.append(gene)

    def amino_acid_frequency(
        self, threshold_unexpected: float = 0.0
    ) -> None:
        """Calculate the frequency of each unique amino acid across all CDS genes in a GeneSet.

        Args:
            threshold_unexpected (float): Percentage of unexpected (not in AA_LIST) amino acids tolerated in a single CDS gene included in the GeneSet (default 0.0 or 0%)
        Populates the following class attributes:
            self.aa_frq (str: float): Frequency of each unique amino acid across all CDS genes in the GeneSet.
        If not populated previously by running amino_acid_counts():
            self.aa_dict (str: int): Counts of each unique amino acid across all CDS genes in the GeneSet.
            self.unexpected_aas (list(str)): list of unexpected amino acids found in the GeneSet.
            self.skipped_unexpected_peptides (List[Gene]): List of CDSGenes in the GeneSet that have more than threshold_unexpected amino acids.
        """
        # Check if GeneSet has any CDS genes
        if len(self.cds_genes) == 0:
            print(f"No Valid CDS genes in {self.id}. Skipping amino acid frequency calculation.")
            return

        print(f"Calculating amino acid frequencies in {self.id}.")
        self.aa_frq: dict[str, float] = {}

        if not hasattr(self, "aa_dict"):
            # If aggregate amino acid counts have not already been calculated, runs amino_acid_counts()
            self.amino_acid_counts(
                threshold_unexpected=threshold_unexpected
            )

        if hasattr(self, "aa_dict") and len(self.skipped_unexpected_peptides)/len(self.cds_genes) < 1:
            total = sum(self.aa_dict.values())
            if total >= 1:
                self.aa_frq = {k: (v / total) for k, v in self.aa_dict.items()}
            elif total < 1:
                print(f"No valid amino acids in {self.id}. Cannot calculate frequencies.")
        else:
            print(f"No valid and precise CDS genes in {self.id} to calculate amino acid frequencies.")

    def RSCU(
        self, threshold_imprecise: float = 0.0
    ) -> None:
        """Calculate the relative synonymous codon usage (RSCU) of each codon across CDS Genes of an entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single CDS gene (default 0.0 or 0%)
        Definitions:
            Synonymous codons: codons that encode the same amino acid
            RSCU_dict: codon count / expected frequency assuming equally abundant synonymous codons
        Populates the following class attributes:
            self.RSCU (str: float): RSCU of each codon across all CDS genes in the GeneSet.
        If not populated previously by running codon_counts() or codon_frequency():
            self.codon_dict (str: int): Counts of each unique codon across all CDS genes in the GeneSet.
            self.imprecise_codons (list(str)): list of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[Gene]): List of CDSGenes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.RSCU_dict: dict[str, float] = dict.fromkeys(CODON_LIST, 0.0)

        if not hasattr(self, "codon_dict"):
            # If aggregate codon counts have not already been calculated, runs codon_counts()
            self.codon_counts(
                threshold_imprecise=threshold_imprecise,
            )

        if hasattr(self, "codon_dict"):
            # create dictionary of the sum of synonymous codon counts for each aminon acid
            expected_counts: dict[str, float] = dict.fromkeys(CODON_LIST, 0)
            for aa in CODON_TABLE.values():
                synonymous_codons = [key for key, value in CODON_TABLE.items() if value == aa]  # list of other codons encoding the same aa
                synonymous_total_count= sum([self.codon_dict[syn_codon] for syn_codon in synonymous_codons])  # total number of synonymous codons present in GeneSet
                expected_counts[aa] = synonymous_total_count / len(synonymous_codons) # expected frequency of codon family members given assumption that all synonymous codons are equally likely to encode the aa

            # calculate the RSCU for each codon
            for codon, count in self.codon_dict.items():
                if count > 0:
                    aa = CODON_TABLE[codon]  # aa encoded by current codon iteration
                    self.RSCU_dict[codon] = self.codon_dict[codon] / expected_counts[aa]

    def tRNA_counts(self) -> None:
        """Calculate the copy numbers of individual tRNA genes by their associated amino acids and (anti)codons.

        Populates the following class attributes:
            self.tRNA_dict_aa (str: int): Counts of tRNA genes by amino acid across all genes in the GeneSet.
            self.tRNA_dict_tcc (str: int): Counts of tRNA genes by their 'tcc' (tRNA complementary codons) across all genes in the GeneSet.
            self.total_tRNA (int): Total number of tRNA genes (not unique) in the GeneSet.
        """
        # Initialize tRNA count dictionaries, skipping stop codons
        self.tRNA_dict_aa: dict[str, int] = {aa: 0 for aa in AA_LIST}
        self.tRNA_dict_tcc: dict[str, int] = {
            tcc: 0 for tcc in CODON_LIST if tcc not in stop_codons
        }

        # Define RegEx pattern for BAKTA-output tRNA gene products
        gene_product_pattern = re.compile(r"tRNA-\w{3}\(\w{3}\)")

        # Populate tRNA count dictionaries by looping through all gene products in GeneSet
        has_tRNAs = False
        for gene in self.genes:
            if gene_product_pattern.match(gene.gene_product):
                has_tRNAs = True
                aa_3 = gene.gene_product.split("-")[1].split("(")[0]
                if aa_3 in AA_CONVERSIONS.keys():
                    aa_1 = AA_CONVERSIONS[aa_3] #new (SeC)
                    self.tRNA_dict_aa[aa_1] += 1
                    anticodon = gene.gene_product.split("(")[1].split(")")[0]
                    tcc = reverse_complement(anticodon)
                    if tcc in self.tRNA_dict_tcc.keys(): #new (UGG)
                        self.tRNA_dict_tcc[tcc] += 1

        # Calculate total tRNA count
        self.total_tRNA: int = sum(self.tRNA_dict_tcc.values())

        if not has_tRNAs:
            print("No tRNA genes found in the GeneSet.")

    def tRNA_frequency(self) -> None:
        """Calculate the frequency of individual tRNA genes (by their associated amino acids and (anti)codons) out of all tRNA genes in the GeneSet.

        Populates the following class attributes:
            self.tRNA_frq_aa (str: int): Frequencies of tRNA genes by amino acid out of all tRNA genes in the GeneSet.
            self.tRNA_frq_tcc (str: int): Frequencies of tRNA genes by their 'tcc' (tRNA complementary codons) out of all tRNA genes in the GeneSet.
        """
        # If tRNA gene counts have not already been calculated, runs tRNA_counts()
        if (
            not hasattr(self, "tRNA_dict_aa")
            or not hasattr(self, "tRNA_dict_tcc")
            or not hasattr(self, "total_tRNA")
        ):
            self.tRNA_counts()

        # Initialize tRNA frequency dictionaries
        self.tRNA_frq_aa: dict[str, float] = dict.fromkeys(
            self.tRNA_dict_aa.keys(), 0.0
        )
        self.tRNA_frq_tcc: dict[str, float] = dict.fromkeys(
            self.tRNA_dict_tcc.keys(), 0.0
        )

        # Calculate frequency of tRNA genes out of total tRNA counts for the GeneSet
        if self.total_tRNA > 0:
            self.tRNA_frq_aa = {
                k: (v / self.total_tRNA) for k, v in self.tRNA_dict_aa.items()
            }
            self.tRNA_frq_tcc = {
                k: (v / self.total_tRNA) for k, v in self.tRNA_dict_tcc.items()
            }


class CodonBiasComparison:
    """Class for calculating codon bias similarity between a virus GeneSet and a host GeneSet.

    Args:
        host_dict (str: int or str: float): Dictionary of codons/amino acids and their counts/frequencies/RSCUs in a host GeneSet.
        virus_dict (str: int or str: float): Dictionary of codons/amino acids and their counts/frequencies/RSCUs in a virus GeneSet.
    """

    def __init__(
        self,
        host_dict: Union[dict[str, int], dict[str, float]],
        virus_dict: Union[dict[str, int], dict[str, float]],
    ) -> None:
        """Initialize class variables and read in an annotated genes file, storing Gene objects and metadata in lists."""
        self.host_dict: Union[dict[str, int], dict[str, float]] = host_dict
        self.host_list: Union[List[int], List[float]] = list(self.host_dict.values())
        self.virus_dict: Union[dict[str, int], dict[str, float]] = virus_dict
        self.virus_list: Union[List[int], List[float]] = list(self.virus_dict.values())

    def linear_regress(self) -> None:
        """Compute linear regression between host and virus codon bias.

        Populates the following class attributes:
            self.slope (float): Slope of the linear regression line between values from input host_dict and virus_dict
            self.intercept (float): Intercept of the linear regression line between values from input host_dict and virus_dict
        """
        lin_regress = np.polyfit(
            np.array(self.host_list), np.array(self.virus_list), 1
        )  # degree is 1 for linear regression
        self.slope: float = float(lin_regress[0])
        self.intercept: float = float(lin_regress[1])

    def calculate_R2(self) -> None:
        """Compute R^2 value between host and virus codon bias using linear regression.

        Populates the following class attributes:
            self.R2 (float): extracts and calculates the R^2 value from linear regression calculation on values from input host_dict and virus_dict
        If not populated previously by running linear_regress():
            self.slope (float): Slope of the linear regression line between values from input host_dict and virus_dict
            self.intercept (float): Intercept of the linear regression line between values from input host_dict and virus_dict
        """
        if not hasattr(self, "slope") and not hasattr(self, "intercept"):
            self.linear_regress()

        if hasattr(self, "slope") and hasattr(self, "intercept"):
            x = np.array(self.host_list)
            y = np.array(self.virus_list)
            y_predicted = self.slope * x + self.intercept
            residuals_sum_sq = np.sum((y - y_predicted) ** 2)
            total_sum_sq = np.sum((y - np.mean(y)) ** 2)
            self.R2: float = float(1 - (residuals_sum_sq / total_sum_sq))

    def cosine_similarity(self):
        """Compute cosine similarity metric between host and virus codon bias.

        Populates the following class attributes:
            self.cos_similarity (float): calculates the cosine similarity between values from input host_dict and virus_dict
        """
        self.cos_similarity: float = float(
            1 - scipy.spatial.distance.cosine(self.host_list, self.virus_list)
        )


class tRNAMetrics:
    """Class for calculating metrics involving tRNA availability.

    Args:
        virus_GeneSet (GeneSet): GeneSet object representing the virus.
        host_GeneSet (GeneSet): GeneSet object representing the host.
    """

    def __init__(self, virus_GeneSet: GeneSet, host_GeneSet: GeneSet) -> None:
        """Initialize class variables."""
        self.virus_GeneSet: GeneSet = virus_GeneSet
        self.host_GeneSet: GeneSet = host_GeneSet
        """Calculate tRNA counts, totals, and frequencies for both virus and host GeneSets."""
        if not hasattr(self.virus_GeneSet, "tRNA_frq_aa"):
            self.virus_GeneSet.tRNA_frequency()
        if not hasattr(self.host_GeneSet, "tRNA_frq_aa"):
            self.host_GeneSet.tRNA_frequency()

    def virus_TAAI(self, include_virus_tRNA: bool = True) -> None:
        """Calculate accordance index between virus amino acid frequency and corresponding tRNA availability. Note that all amino acids are included in the correlation.

        Args:
            include_virus_tRNA (bool): Whether to additionally calculate an accordance metric that accounts for tRNA gene counts from virus, in addition to that of host (default is True).

        Populates the following class attributes:
            self.virusTAAI_hosttRNA (float): Spearman rank correlation coefficient between host tRNA gene copy frequencies and corresponding viral amino acid frequencies.
            self.virusTAAI_totaltRNA (float): Attribute created and populated only if include_virus_tRNA argument is set to True. Spearman rank correlation coefficient between total tRNA gene copy frequencies (virus and host) and corresponding viral amino acid frequencies.
        """
        # Generate amino acid frequencies for virus GeneSet if not already existent
        if not hasattr(self.virus_GeneSet, "aa_frq"):
            self.virus_GeneSet.amino_acid_frequency()

        # Perform Spearman Rank correlation between virus amino acid frequency and host tRNA availability
        sorted_keys = sorted(self.virus_GeneSet.aa_frq)
        virus_aa_frq_values = [self.virus_GeneSet.aa_frq[key] for key in sorted_keys]
        host_tRNA_frq_aa_values = [
            self.host_GeneSet.tRNA_frq_aa[key] for key in sorted_keys
        ]
        res = scipy.stats.spearmanr(virus_aa_frq_values, host_tRNA_frq_aa_values)
        self.virusTAAI_hosttRNA: float = res.statistic

        # If specified, perform Spearman Rank correlation between virus amino acid frequency and total tRNA availability
        if include_virus_tRNA is True:
            total_tRNA_dict_aa = {
                key: self.host_GeneSet.tRNA_dict_aa.get(key, 0)
                + self.virus_GeneSet.tRNA_dict_aa.get(key, 0)
                for key in set(self.host_GeneSet.tRNA_dict_aa)
                | set(self.virus_GeneSet.tRNA_dict_aa)
            }
            total_virocell_tRNA = sum(total_tRNA_dict_aa.values())
            total_tRNA_frq_aa = {
                k: (v / total_virocell_tRNA) for k, v in total_tRNA_dict_aa.items()
            }
            total_tRNA_frq_aa_values = [total_tRNA_frq_aa[key] for key in sorted_keys]
            res = scipy.stats.spearmanr(virus_aa_frq_values, total_tRNA_frq_aa_values)
            self.virusTAAI_totaltRNA: float = res.statistic

    def virus_TCAI(
        self, skip_nondeg_codons: bool = True, include_virus_tRNA: bool = True
    ) -> None:
        """Calculate accordance index between virus codon frequency and corresponding tRNA availability.

        Args:
            skip_nondeg_codons (bool): Whether to omit non-degenerate codons (codons whose encoded amino acid is specific to one codon alone) from the accordance calculation (default is True). Note stop codons are inherently skipped because they have no associated tRNA.
            include_virus_tRNA (bool): Whether to additionally calculate an accordance metric that accounts for tRNA gene counts from virus, in addition to that of host (default is True).

        Populates the following class attributes:
            self.virusTCAI_hosttRNA (float): Spearman rank correlation coefficient between host tRNA gene copy frequencies and corresponding viral codon frequencies.
            self.virusTCAI_totaltRNA (float): Attribute created and populated only if include_virus_tRNA argument is set to True. Spearman rank correlation coefficient between total tRNA gene copy frequencies (virus and host) and corresponding viral codon frequencies.
        """
        # Generate codon frequencies for virus GeneSet if not already existent
        if not hasattr(self.virus_GeneSet, "codon_frq"):
            self.virus_GeneSet.codon_frequency()

        # Define tRNA dictionaries, skipping non-degenerate codons and including virus tRNAs as specified
        if skip_nondeg_codons is True:
            # Sort and remove non-degenerate and stop codon keys
            irrelevant_codons = non_degenerate_codons + stop_codons
            sorted_keys = [
                key
                for key in sorted(self.virus_GeneSet.codon_frq)
                if key not in irrelevant_codons
            ]
            # Prepare codon and tRNA lists for spearman rank
            virus_codon_frq_values = [
                self.virus_GeneSet.codon_frq[key] for key in sorted_keys
            ]
            host_tRNA_frq_tcc_values = [
                self.host_GeneSet.tRNA_frq_tcc[key] for key in sorted_keys
            ]
            # Prepare tRNA dictionaries for total tRNA accordance metric if specified
            virus_tRNA_dict_tcc = (
                {
                    k: v
                    for k, v in self.virus_GeneSet.tRNA_dict_tcc.items()
                    if k not in irrelevant_codons
                }
                if include_virus_tRNA
                else None
            )
            host_tRNA_dict_tcc = (
                {
                    k: v
                    for k, v in self.host_GeneSet.tRNA_dict_tcc.items()
                    if k not in irrelevant_codons
                }
                if include_virus_tRNA
                else None
            )
        else:
            # Sort and remove stop codon keys only
            sorted_keys = [
                key
                for key in sorted(self.virus_GeneSet.codon_frq)
                if key not in stop_codons
            ]
            # Prepare codon and tRNA lists for spearman rank
            virus_codon_frq_values = [
                self.virus_GeneSet.codon_frq[key] for key in sorted_keys
            ]
            host_tRNA_frq_tcc_values = [
                self.host_GeneSet.tRNA_frq_tcc[key] for key in sorted_keys
            ]
            # Prepare tRNA dictionaries for total tRNA accordance metric if specified
            virus_tRNA_dict_tcc = (
                {
                    k: v
                    for k, v in self.virus_GeneSet.tRNA_dict_tcc.items()
                    if k not in stop_codons
                }
                if include_virus_tRNA
                else None
            )
            host_tRNA_dict_tcc = (
                {
                    k: v
                    for k, v in self.host_GeneSet.tRNA_dict_tcc.items()
                    if k not in stop_codons
                }
                if include_virus_tRNA
                else None
            )

        # Perform Spearman Rank correlation between virus codon frequency and host tRNA availability
        res = scipy.stats.spearmanr(virus_codon_frq_values, host_tRNA_frq_tcc_values)
        self.virusTCAI_hosttRNA: float = res.statistic

        # If specified, perform Spearman Rank correlation between virus codon frequency and total tRNA availability
        if (
            include_virus_tRNA is True
            and host_tRNA_dict_tcc is not None
            and virus_tRNA_dict_tcc is not None
        ):
            total_tRNA_dict_tcc = {
                key: host_tRNA_dict_tcc.get(key, 0) + virus_tRNA_dict_tcc.get(key, 0)
                for key in set(host_tRNA_dict_tcc) | set(virus_tRNA_dict_tcc)
            }
            total_virocell_tRNA = sum(total_tRNA_dict_tcc.values())
            total_tRNA_frq_tcc_values = [
                total_tRNA_dict_tcc[key] / total_virocell_tRNA for key in sorted_keys
            ]
            res = scipy.stats.spearmanr(
                virus_codon_frq_values, total_tRNA_frq_tcc_values
            )
            self.virusTCAI_totaltRNA: float = res.statistic

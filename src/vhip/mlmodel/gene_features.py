"""Contains classes to compute gene-level features.

This module provides:
- Gene: calculate codon counts, amino acid counts, and imprecise codon counts for a single gene
- GeneSet: calculate codon counts, codon frequency, RSCU (relative synonymous codon usage), amino acid counts, and amino acid frequency for a set of genes in an annotated file
- CodonBiasComparison: compare codon bias measurements between two gene sets using linear regression (slope, R^2) and cosine similarity
"""

import os
from typing import List, Union

import scipy  # pyright: ignore[reportMissingTypeStubs]

from .read_sequence import read_annotated_genes

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
AA_LIST = list(CODON_TABLE.values())
stop_codons = [codon for codon, aa in CODON_TABLE.items() if aa == "_"]
non_degenerate_codons = [codon for codon, aa in CODON_TABLE.items() if list(CODON_TABLE.values()).count(aa) == 1]

# Amino acid abbreviations conversions
AA_CONVERSIONS = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V'
}


# Define Gene class
class Gene:
    """Class representing a gene.

    Args:
        gene_seq (str): The nucleotide sequence of the gene.
        codon_length (int): Length of 1 codon (default is 3).
        gene_id (str): The gene ID from annotations (default is an empty string).
        gene_product (str): The product of the gene from annotations (default is an empty string).
    """

    def __init__(
        self,
        gene_seq: str,
        codon_length: int = 3,
        gene_id: str = "",
        gene_product: str = "",
    ) -> None:
        """Initialize class variables."""
        if len(gene_seq) % codon_length == 0:
            self.seq = gene_seq
            self.codon_length = codon_length
            self.n_codons = len(gene_seq) / codon_length
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
        self.number_imprecise_codons: int = 0
        self.codon_dict = dict.fromkeys(CODON_LIST, 0)

        for i in range(0, len(self.seq), self.codon_length):
            codon = self.seq[i : i + self.codon_length]
            if codon in self.codon_dict.keys():
                self.codon_dict[codon] += 1
            else:
                self.number_imprecise_codons += 1

        self.percent_imprecise_codons: float = (
            self.number_imprecise_codons / self.n_codons
        )

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

        for i in range(0, len(self.seq), self.codon_length):
            codon = self.seq[i : i + self.codon_length]
            for j in range(self.codon_length):
                if codon[j] == "G" or codon[j] == "C":
                    if j == 0:
                        gc1 += 1
                    elif j == 1:
                        gc2 += 1
                    elif j == 2:
                        gc3 += 1

        self.GC1 = gc1 / self.n_codons
        self.GC2 = gc2 / self.n_codons
        self.GC3 = gc3 / self.n_codons


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
        self.id = gene_file
        self.genes: List[Gene] = []
        self.skipped_genes: List[str] = []
        readout = read_annotated_genes(gene_file)

        for out in range(len(readout[0])):
            try:
                self.genes.append(
                    Gene(
                        gene_seq=str(readout[0][out]),
                        gene_id=str(readout[1][out]),
                        gene_product=str(readout[2][out]),
                    )
                )
            except Exception:
                self.skipped_genes.append(str(readout[1][out]))
        percent_skipped = len(self.skipped_genes) / len(readout[0]) * 100
        print(
            f"{percent_skipped}% ({len(self.skipped_genes)}/{len(readout[0])}) of genes skipped. Expect on average ~2% and ~3% of virus and host genes (respectively) to be skipped on the basis of non-divisibility by codon length."
        )

    def codon_counts(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Aggregate the counts for each unique codon and imprecise codons across an entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0%)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50%)
        Populates the following class attributes:
            self.codon_dict (str: int): Counts of each unique codon across all genes in the GeneSet.
            self.imprecise_codons (int): Total number of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List(str)): IDs of genes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.codon_dict: dict[str, int] = dict.fromkeys(CODON_LIST, 0)
        self.imprecise_codons: int = 0
        self.skipped_imprecise_genes: List[str] = []

        counter = 0
        for gene in self.genes:
            counter += 1
            print(f"Analyzing gene {counter} of {len(self.genes)}")
            gene.calculate_codon_counts()
            self.imprecise_codons += gene.number_imprecise_codons
            if gene.percent_imprecise_codons <= threshold_imprecise:
                for key, val in gene.codon_dict.items():
                    self.codon_dict[key] += val
            else:
                self.skipped_imprecise_genes.append(gene.gene_id)

        if (
            len(self.skipped_imprecise_genes) / len(self.genes)
            > threshold_skipped_genes
        ):
            raise Exception(
                f"Too many skipped genes. {len(self.skipped_imprecise_genes)} genes have > {threshold_imprecise} imprecise codons."
            )
        elif len(self.skipped_imprecise_genes) > 0:
            print(
                f"Skipped {len(self.skipped_imprecise_genes)} genes with too many imprecise codons"
            )

    def codon_frequency(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Calculate the frequency of each unique codon in an entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0%)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50%)
        Populates the following class attributes:
            self.codon_frq (str: float): Frequency of each unique codon across all genes in the GeneSet.
        If not populated previously by running codon_counts():
            self.codon_dict (str: int): Counts of each unique codon across all genes in the GeneSet.
            self.imprecise_codons (int): Total number of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[str]): IDs of genes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.codon_frq: dict[str, float] = {}

        if not hasattr(self, "codon_dict"):
            # If aggregate codon counts have not already been calculated, runs codon_counts()
            self.codon_counts(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

        if hasattr(self, "codon_dict"):
            total = sum(self.codon_dict.values())
            self.codon_frq = {k: (v / total) for k, v in self.codon_dict.items()}

    def amino_acid_counts(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Calculate the counts of each unique amino acid encoded by an entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0%)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50%)
        Populates the following class attributes:
            self.aa_dict (str: int): Counts of each unique amino acid across all genes in the GeneSet.
        If not populated previously by running codon_counts():
            self.codon_dict (str: int): Counts of each unique codon across all genes in the GeneSet.
            self.imprecise_codons (int): Total number of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[str]): IDs of genes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.aa_dict: dict[str, int] = dict.fromkeys(AA_LIST, 0)

        if not hasattr(self, "codon_dict"):
            # If aggregate codon counts have not already been calculated, runs codon_counts()
            self.codon_counts(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

        if hasattr(self, "codon_dict"):
            for codon in self.codon_dict.keys():
                aa = CODON_TABLE[codon]
                self.aa_dict[aa] += self.codon_dict[codon]

    def amino_acid_frequency(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Calculate the frequency of each unique amino acid encoded by an entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0%)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50%)
        Populates the following class attributes:
            self.aa_frq (str: float): Frequency of each unique amino acid across all genes in the GeneSet.
        If not populated previously by running amino_acid_counts() and codon_counts():
            self.aa_dict (str: int): Counts of each unique amino acid across all genes in the GeneSet.
            self.codon_dict (str: int): Counts of each unique codon across all genes in the GeneSet.
            self.imprecise_codons (int): Total number of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[str]): IDs of genes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.aa_frq: dict[str, float] = {}

        if not hasattr(self, "aa_dict"):
            # If aggregate amino acid counts have not already been calculated, runs amino_acid_counts()
            self.amino_acid_counts(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

        if hasattr(self, "aa_dict"):
            total = sum(self.aa_dict.values())
            self.aa_frq = {k: (v / total) for k, v in self.aa_dict.items()}

    def RSCU(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Calculate the relative synonymous codon usage (RSCU) of each codon in entire GeneSet.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0%)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50%)
        Definitions:
            Synonymous codons: codons that encode the same amino acid
            RSCU_dict: codon count / expected frequency (given assumption of equally used synonymous codons)
        Populates the following class attributes:
            self.RSCU (str: float): RSCU of each codon across all genes in the GeneSet.
        If not populated previously by running codon_counts() or codon_frequency():
            self.codon_dict (str: int): Counts of each unique codon across all genes in the GeneSet.
            self.imprecise_codons (int): Total number of imprecise codons found in the GeneSet.
            self.skipped_imprecise_genes (List[str]): IDs of genes in the GeneSet that have more than threshold_imprecise codons.
        """
        self.RSCU_dict: dict[str, float] = dict.fromkeys(CODON_LIST, 0.0)

        if not hasattr(self, "codon_dict"):
            # If aggregate codon counts have not already been calculated, runs codon_counts()
            self.codon_counts(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

        if hasattr(self, "codon_dict"):
            for codon, count in self.codon_dict.items():
                if count != 0.0:
                    aa = CODON_TABLE[codon]  # aa encoded by current codon iteration
                    synonymous_codons = [
                        key for key, value in CODON_TABLE.items() if value == aa
                    ]  # list of other codons encoding the same aa
                    synonymous_total_count = sum(
                        [self.codon_dict[syn_codon] for syn_codon in synonymous_codons]
                    )  # total number of synonymous codons present in GeneSet
                    expected_frequency = (
                        synonymous_total_count / len(synonymous_codons)
                    )  # expected frequency of current codon given all assumption all synonymous codons are equally likely to encode the aa
                    self.RSCU_dict[codon] = self.codon_dict[codon] / expected_frequency

    def tRNA_counts(self, skip_nondeg_codons: bool = True) -> None:
        """Calculate the copy numbers of individual tRNA genes by their associated amino acids and (anti)codons.

        Args:
            skip_nondeg_codons (bool): Whether to include non-degenerate codons (codons whose encoded amino acid is specific to one codon alone) in the list of codons to skip (default is True). Note stop codons are automatically skipped because they have no associated amino acid/tRNA.

        Populates the following class attributes:
            self.tRNA_dict_aa (str: int): Counts of tRNA genes by amino acid across all genes in the GeneSet.
            self.tRNA_dict_tcc (str: int): Counts of tRNA genes by their 'tcc' (tRNA complementary codons) across all genes in the GeneSet.
        """
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
            self.lin_regress (LinregresResult class from scipy.stats._stats_mstats_common): output of running linear regression between the values from input host_dict and virus_dict
        """
        self.lin_regress: scipy.stats._stats_mstats_common.LinregressResult = (
            scipy.stats.linregress(self.host_list, self.virus_list)
        )

    def calculate_slope(self) -> None:
        """Compute slope between host and virus codon bias using linear regression.

        Populates the following class attributes:
            self.slope (float): extracts slope from linear regression calculation on values from input host_dict and virus_dict
        If not populated previously by running linear_regress():
            self.lin_regress (LinregresResult class from scipy.stats._stats_mstats_common): output of running linear regression between the values from input host_dict and virus_dict
        """
        if not hasattr(self, "lin_regress"):
            self.linear_regress()

        if hasattr(self, "lin_regress"):
            self.slope: float = float(
                self.lin_regress.slope
            )  # the first value from the output of scipy.stats.linregress is slope

    def calculate_R2(self) -> None:
        """Compute R^2 value between host and virus codon bias using linear regression.

        Populates the following class attributes:
            self.R2 (float): extracts and calculates the R^2 value from linear regression calculation on values from input host_dict and virus_dict
        If not populated previously by running linear_regress():
            self.lin_regress (LinregresResult class from scipy.stats._stats_mstats_common): output of running linear regression between the values from input host_dict and virus_dict
        """
        if not hasattr(self, "lin_regress"):
            self.linear_regress()

        if hasattr(self, "lin_regress"):
            self.R2: float = float(
                self.lin_regress.rvalue**2
            )  # the third value from the output of scipy.stats.linregress is the r-value

    def cosine_similarity(self):
        """Compute cosine similarity metric between host and virus codon bias.

        Populates the following class attributes:
            self.cos_similarity (float): calculates the cosine similarity between values from input host_dict and virus_dict
        """
        self.cos_similarity: float = float(
            1 - scipy.spatial.distance.cosine(self.host_list, self.virus_list)
        )

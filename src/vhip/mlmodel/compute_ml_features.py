"""Methods and functions to compute signals of coevolutions between viruses and putative hosts.

Those signals are necessary for virus-host predictions.
"""

import itertools
import os
from multiprocessing import Pool
from typing import List, Optional

import pandas as pd  # pyright: ignore[reportMissingTypeStubs]

from .gene_features import CodonBiasComparison, GeneSet
from .genomes_features import HomologyMatch, KmerProfile, d2Distance
from .read_sequence import read_headers, read_sequence


class Pairs:
    """Class to store coevolution signal values for a virus-host pair.

    Args:
        virus (str): filename of virus of interest
        host (str): filename of host of interest
    """

    def __init__(self, virus: str, host: str) -> None:
        """Initialize class variables."""
        self.virus = virus
        self.host = host

        # Genome-level features
        self.GCdifference: float
        self.k3dist: float
        self.k6dist: float
        self.homology_hit: bool

        # Gene-level features
        self.codons_comparison: CodonBiasComparison
        self.aa_comparison: CodonBiasComparison
        self.RSCU_comparison: CodonBiasComparison

        self.interaction: int


class ComputeFeatures:
    """Class organizing all methods to compute all virus-host coevolution signals.

    Args:
        virus_genome_dir (str): Path to the directory containing viruses genome fasta files. Each file should contain an unique virus.
        host_genome_dir (str): Path to the directory containing host genome fasta files.  Each file should contain an unique host species/OTUs.
        virus_gene_dir (str): Path to the directory containing viruses gene fasta files. Each file should contain annotated genes for a unique virus.
        host_gene_dir (str): Path to the directory containing host gene fasta files.  Each file should contain annotated genes for a unique host species/OTUs.
        genome_ext (str): Extension used for genome fasta files. Default is "fasta".
        gene_ext (str): Extension used for gene fasta files. Default is "ffn".
        pairs_of_interest (str): Pathway to file containing virus-host pairs of interest. Optional.
    """

    def __init__(
        self,
        virus_genome_dir: str,
        host_genome_dir: str,
        virus_gene_dir: str,
        host_gene_dir: str,
        genome_ext: str = "fasta",
        gene_ext: str = "ffn",
        pairs_of_interest: Optional[str] = None,

    ) -> None:
        """Initialize class variables."""
        self.virus_genome_dir = virus_genome_dir
        self.host_genome_dir = host_genome_dir
        self.genome_ext = genome_ext
        self.virus_gene_dir = virus_gene_dir
        self.host_gene_dir = host_gene_dir
        self.gene_ext = gene_ext
        self.pairs_of_interest = pairs_of_interest

        self.features_df = None

    def add_blastn_files(self, blastn_path: str, spacer_path: str):
        """Add blastn path for viruses against hosts and viruses against spacers.

        Args:
            blastn_path (str): Path to blastn output of viruses against hosts
            spacer_path (str): Path to blastn output of viruses against spacers
        """
        self.blastn_path = blastn_path
        self.spacer_path = spacer_path

    def do_setup(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ):
        """Calls other methods to setup.

        The setup process includes determining all possible virus-host pairs, get fasta headers, read and process blastn_output, compute GC content and k-mer profiles, and generate dictionaries of codon, amino acid, and synonymous codon usage frequencies.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0% - see paper methods for threshold default determination)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50% - see paper methods for threshold default determination)
        """
        print("SETUP - ...indexing genome fasta filenames for viruses and hosts...")
        self.list_genome_files()

        print(
            "SETUP - ...indexing annotated gene fasta filenames for viruses and hosts..."
        )
        self.list_gene_files()

        print("SETUP - ...initialize all pairs...")
        if self.pairs_of_interest:
            self.determine_custom_pairs(self.pairs_of_interest)
        else:
            self.determine_pairs()

        print("SETUP - ...getting fasta headers...")
        self.get_headers()

        print("SETUP - ...process blastn and spacers output...")
        self.process_blastn()
        self.process_spacers()
        self.homology_matches = HomologyMatch(self.blastn, self.spacers)

        print("SETUP - ...calculate GC content and k-mer profiles...")
        self.generate_kmer_profiles()

        print("SETUP - ...calculate codon and amino acid profiles...")
        self.generate_codon_frq(threshold_imprecise, threshold_skipped_genes)
        self.generate_aa_frq(threshold_imprecise, threshold_skipped_genes)
        self.generate_RSCU(threshold_imprecise, threshold_skipped_genes)

    def list_genome_files(self):
        """List all genome fasta file in the virus and host genome directories."""
        self.virus_genome_filenames = [
            f
            for f in os.listdir(self.virus_genome_dir)
            if f.endswith("." + self.genome_ext)
        ]
        self.host_genome_filenames = [
            f
            for f in os.listdir(self.host_genome_dir)
            if f.endswith("." + self.genome_ext)
        ]
        self.all_genome_files = self.virus_genome_filenames + self.host_genome_filenames

    def list_gene_files(self):
        """List all gene fasta file in the virus and host annotated gene file directories."""
        self.virus_gene_filenames = [
            f
            for f in os.listdir(self.virus_gene_dir)
            if f.endswith("." + self.gene_ext)
        ]
        self.host_gene_filenames = [
            f for f in os.listdir(self.host_gene_dir) if f.endswith("." + self.gene_ext)
        ]
        self.all_gene_files = self.virus_gene_filenames + self.host_gene_filenames

    def determine_pairs(self):
        """Determine all possible virus-host pairs.

        This assume that each virus should be tested against each host.
        """
        total_interactions = len(self.virus_genome_filenames) * len(
            self.host_genome_filenames
        )
        print(f"-------> There are {len(self.virus_genome_filenames)} viral sequences")
        print(f"-------> There are {len(self.host_genome_filenames)} host sequences")
        print(f"-------> Total number of interactions: {total_interactions}")

        # determine all virus-host pair possible (every host is going to be considered for every virus of interest)
        virus_inter: List[str] = list(
            itertools.chain.from_iterable(
                itertools.repeat(x, len(self.host_genome_filenames))
                for x in self.virus_genome_filenames
            )
        )
        host_inter = self.host_genome_filenames * len(self.virus_genome_filenames)

        # create list of Pairs
        self.pairs = []
        for pair in list(zip(virus_inter, host_inter)):
            virus = pair[0]
            host = pair[1]
            self.pairs.append(Pairs(virus, host))

    def determine_custom_pairs(self, custom_pairs: str):
        """The input pair file needs to be virus first then host. Must be separated by tabs."""
        self.pairs: List[Pairs] = []
        print("reading pairs file")

        with open(custom_pairs, "r") as f:
            lines = [line.rstrip() for line in f]
            for pair in lines:
                split = pair.split("\t")
                virus = split[0]
                host = split[1]

                if (
                    virus in self.virus_genome_filenames
                    and host in self.host_genome_filenames
                ):
                    self.pairs.append(Pairs(virus, host))
                else:
                    print(
                        "Warning: file pair "
                        + virus
                        + " and "
                        + host
                        + " are not all present. Please double check."
                    )

    def get_headers(self):
        """Retrieve headers from genome fasta files.

        This is used to map where each blast hit came from between viruses and hosts.
        """
        header_filename: dict[str, str] = {}

        for virus in self.virus_genome_filenames:
            path = self.virus_genome_dir + virus
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = virus

        for host in self.host_genome_filenames:
            path = self.host_genome_dir + host
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = host

        self.headers = header_filename

    def process_blastn(self):
        """Process blastn file of viruses against hosts and store data in a dictionary.

        Keys of the dictionary are viruses. For each key, there is a list of hosts that had a homology hit with that virus.
        """
        self.blastn: dict[str, list[str]] = {}

        with open(self.blastn_path, "r") as f:
            lines = [line.rstrip() for line in f]
            for line in lines:
                split = line.split("\t")
                # get accessions (blastn naming system)
                virus_acc = split[0]
                host_acc = split[1]
                # switch to naming by filenames
                virus = self.headers.get(virus_acc, "NA")
                host = self.headers.get(host_acc, "NA")
                # add virus host relation to dictionary
                if virus not in self.blastn:
                    self.blastn[virus] = [host]
                elif host not in self.blastn[virus]:
                    self.blastn[virus].append(host)

    def process_spacers(self):
        """Process blastn file of viruses against spacers and store data in a dictionary.

        Key of the dictionary are viruses. For each key, there is a list of hosts that had a homology hit with that virus.
        """
        self.spacers: dict[str, list[str]] = {}

        with open(self.spacer_path, "r") as f:
            lines = [line.rstrip() for line in f]
            for line in lines:
                split = line.split("\t")
                # get filename for virus
                virus_acc = split[0]
                virus = self.headers[virus_acc]
                # get filename for host
                tmp = split[1].split("_")
                host_partial = tmp[0] + "_" + tmp[1]
                host = list(
                    filter(
                        lambda x: x.startswith(host_partial), self.host_genome_filenames
                    )
                )
                if virus not in self.spacers:
                    self.spacers[virus] = host
                else:
                    self.spacers[virus].append(host[0])

    def generate_kmer_profiles(self):
        """Generate kmer profiles for viruses and hosts, and compute GC content.

        This will be done for each fasta files in the virus and host directories.
        """
        self.GCcontent = dict.fromkeys(self.all_genome_files)
        self.k3profiles = dict.fromkeys(self.all_genome_files)
        self.k6profiles = dict.fromkeys(self.all_genome_files)

        for virus in self.virus_genome_filenames:
            path = self.virus_genome_dir + virus
            seq = read_sequence(path)

            seq_profile = KmerProfile(seq, k=1)
            seq_profile.generate_profile()
            self.GCcontent[virus] = seq_profile.GCcontent

            seq_profile = KmerProfile(seq, k=3)
            seq_profile.generate_profile()
            self.k3profiles[virus] = seq_profile

            seq_profile = KmerProfile(seq, k=6)
            seq_profile.generate_profile()
            self.k6profiles[virus] = seq_profile

        for host in self.host_genome_filenames:
            path = self.host_genome_dir + host
            seq = read_sequence(path)

            seq_profile = KmerProfile(seq, k=1)
            seq_profile.generate_profile()
            self.GCcontent[host] = seq_profile.GCcontent

            seq_profile = KmerProfile(seq, k=3)
            seq_profile.generate_profile()
            self.k3profiles[host] = seq_profile

            seq_profile = KmerProfile(seq, k=6)
            seq_profile.generate_profile()
            self.k6profiles[host] = seq_profile

    def generate_codon_aa_counts(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Set up GeneSet objects for each virus and host gene files, populate their codon_dict and aa_dict (counts) variables, and store those dictionaries in ComputeFeatures variables."""
        self.virus_GeneSets = {
            virus: GeneSet(os.path.join(self.virus_gene_dir, virus))
            for virus in self.virus_gene_filenames
        }
        self.host_GeneSets = {
            host: GeneSet(os.path.join(self.host_gene_dir, host))
            for host in self.host_gene_filenames
        }
        self.codon_counts = dict.fromkeys(self.all_gene_files)
        self.aa_counts = dict.fromkeys(self.all_gene_files)

        for virus, virus_GeneSet in self.virus_GeneSets.items():
            virus_GeneSet.amino_acid_counts(
                threshold_imprecise, threshold_skipped_genes
            )
            self.codon_counts[virus] = virus_GeneSet.codon_dict
            self.aa_counts[virus] = virus_GeneSet.aa_dict

        for host, host_GeneSet in self.host_GeneSets.items():
            host_GeneSet.amino_acid_counts(threshold_imprecise, threshold_skipped_genes)
            self.codon_counts[host] = host_GeneSet.codon_dict
            self.aa_counts[host] = host_GeneSet.aa_dict

    def generate_codon_frq(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Generate profile of the frequency of each unique codon in every virus and host GeneSet.

        This will be compiled from all genes in each .ffn file in the virus and host genes files directories.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0% - see paper methods for threshold default determination)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50% - see paper methods for threshold default determination)
        """
        self.codon_frqs: dict[str, dict[str, float]] = {
            key: {} for key in self.all_gene_files
        }

        if not hasattr(self, "codon_counts"):
            # If aggregate codon counts have not already been calculated for the GeneSets, runs generate_codon_aa_counts()
            self.generate_codon_aa_counts(threshold_imprecise, threshold_skipped_genes)

        for virus, virus_GeneSet in self.virus_GeneSets.items():
            virus_GeneSet.codon_frequency(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.codon_frqs[virus] = virus_GeneSet.codon_frq

        for host, host_GeneSet in self.host_GeneSets.items():
            host_GeneSet.codon_frequency(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.codon_frqs[host] = host_GeneSet.codon_frq

    def generate_aa_frq(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Generate profile of the frequency of each unique amino acid in every virus and host GeneSet.

        This will be compiled from all genes in each .ffn file in the virus and host genes files directories.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0% - see paper methods for threshold default determination)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50% - see paper methods for threshold default determination)
        """
        self.aa_frqs: dict[str, dict[str, float]] = {
            key: {} for key in self.all_gene_files
        }

        if not hasattr(self, "aa_counts"):
            # If aggregate amino acid counts have not already been calculated for the GeneSets, runs generate_codon_aa_counts()
            self.generate_codon_aa_counts(threshold_imprecise, threshold_skipped_genes)

        for virus, virus_GeneSet in self.virus_GeneSets.items():
            virus_GeneSet.amino_acid_frequency(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.aa_frqs[virus] = virus_GeneSet.aa_frq

        for host, host_GeneSet in self.host_GeneSets.items():
            host_GeneSet.amino_acid_frequency(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.aa_frqs[host] = host_GeneSet.aa_frq

    def generate_RSCU(
        self, threshold_imprecise: float = 0.0, threshold_skipped_genes: float = 0.5
    ) -> None:
        """Generate profile of the relative synonymous codon usage (RSCU) of each codon in every virus and host GeneSet.

        This will be compiled from all genes in each .ffn file in the virus and host genes files directories.

        Args:
            threshold_imprecise (float): Percentage of imprecise (non-ATGC) codons tolerated in a single gene (default 0.0 or 0% - see paper methods for threshold default determination)
            threshold_skipped_genes (float): Tolerated percentage of valid (codon length divisible) genes in GeneSet that have more than threshold_imprecise codons (default 0.5 or 50% - see paper methods for threshold default determination)
        """
        self.RSCU: dict[str, dict[str, float]] = {
            key: {} for key in self.all_gene_files
        }

        if not hasattr(self, "codon_counts"):
            # If aggregate codon counts have not already been calculated for the GeneSets, runs generate_codon_aa_counts()
            self.generate_codon_aa_counts(threshold_imprecise, threshold_skipped_genes)

        for virus, virus_GeneSet in self.virus_GeneSets.items():
            virus_GeneSet.RSCU(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.RSCU[virus] = virus_GeneSet.RSCU_dict

        for host, host_GeneSet in self.host_GeneSets.items():
            host_GeneSet.RSCU(
                threshold_imprecise=threshold_imprecise,
                threshold_skipped_genes=threshold_skipped_genes,
            )

            self.RSCU[host] = host_GeneSet.RSCU_dict

    def run_parallel(self, num_procs: int = 6):
        """Run multiple process of the compute_feature method.

        This significantly improves run time by using multiple CPU cores.

        Args:
            num_procs (int): Number of core to be used.
        """
        with Pool(num_procs) as pool:
            results = pool.map(self.compute_feature, self.pairs)

        self.computed_pairs = results

    def compute_feature(self, pair: Pairs) -> Pairs:
        """Compute all virus-host coevolution signals needed to predict interaction.

        Args:
            pair (Dataclass): virus-host Pairs dataclass.
        """
        print(f"-------> current pair: {pair.virus} | {pair.host}")
        pair.homology_hit = self.homology_matches.match(pair.virus, pair.host)

        k3distance = d2Distance(self.k3profiles[pair.virus], self.k3profiles[pair.host])  # pyright: ignore
        k3distance.distance()
        pair.k3dist = k3distance.dist  # pyright: ignore

        k6distance = d2Distance(self.k6profiles[pair.virus], self.k6profiles[pair.host])  # pyright: ignore
        k6distance.distance()
        pair.k6dist = k6distance.dist  # pyright: ignore

        pair.GCdifference = self.GCcontent[pair.virus] - self.GCcontent[pair.host]  # pyright: ignore

        # Create CodonBiasComparison objects for the pair
        virus = pair.virus.replace(
            self.genome_ext, self.gene_ext
        )  # first get the correct key based on genes file name (virus)
        host = pair.host.replace(
            self.genome_ext, self.gene_ext
        )  # first get the correct key based on genes file name (host)
        pair.codons_comparison = CodonBiasComparison(
            self.codon_frqs[virus],
            self.codon_frqs[host],
        )
        pair.aa_comparison = CodonBiasComparison(
            self.aa_frqs[virus],
            self.aa_frqs[host],
        )
        pair.RSCU_comparison = CodonBiasComparison(
            self.RSCU[virus],
            self.RSCU[host],
        )
        # for each of the above CodonBiasComparison objects, compute all comparisons (R2, slope, cosine similarity)
        for comparison in [
            pair.codons_comparison,
            pair.aa_comparison,
            pair.RSCU_comparison,
        ]:
            comparison.calculate_slope()
            comparison.calculate_R2()
            comparison.cosine_similarity()

        return pair

    def convert_to_dataframe(self):
        """Convert the list of Pairs into a Pandas dataframe."""
        pairs: List[str] = []

        k3dist: List[float] = []
        k6dist: List[float] = []
        GCdiff: List[float] = []
        Homology: List[int] = []

        for pair in self.computed_pairs:
            virus_host = str(pair.virus + ":" + pair.host)
            pairs.append(virus_host)

            k3dist.append(pair.k3dist)
            k6dist.append(pair.k6dist)
            GCdiff.append(pair.GCdifference)
            Homology.append(int(pair.homology_hit))

        self.features_df = pd.DataFrame(
            list(zip(pairs, GCdiff, k3dist, k6dist, Homology)),
            columns=["pairs", "GCdiff", "k3dist", "k6dist", "Homology"],
        )
        self.features_df = self.features_df.set_index("pairs")  # pyright: ignore[reportUnknownMemberType]

    def save_features(self, filename: str):
        """Save computed features as a tsv file.

        Args:
            filename (str): Name to be used when exporting features as a tsv file.
        """
        if self.features_df is None:
            self.convert_to_dataframe()

        self.features_df.to_csv(filename, sep="\t")  # pyright: ignore

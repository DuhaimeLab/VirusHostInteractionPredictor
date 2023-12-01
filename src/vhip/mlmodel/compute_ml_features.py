"""This script calls methods and functions to compute signals of coevolutions.
Those signals are necessary for virus0-host predictions.
"""

from dataclasses import dataclass
import itertools
import os
import pandas as pd

from multiprocessing import Pool

from .genomes_features import KmerProfile, d2Distance, HomologyMatch
from .read_sequence import read_sequence, read_headers


@dataclass
class Pairs:
    """Dataclass to store coevolution signals values for a specific virus-host pair.
    To initiaize this class, only the filename of the virus and host are needed.
    This is needed for the ComputeFeatures class below.

    :param virus: Filename of virus of interest
    :type virus: str
    :param host: Filename of host of interest
    :type host: str
    """

    virus: str
    host: str
    interaction: int = None

    # genome level features
    GCdifference: float = None
    k3dist: float = None
    k6dist: float = None

    homology_hit: bool = False

    # TODO: gene level features


class ComputeFeatures:
    """The ComputeFeatures class takes for input the virus directory (location of all viral fasta)
    and the host directory (location of all hosts fasta). Its purpose is to determine all possible
    pairs to be computed, and calls on other classes to calculate signals of co-evolution.
    In other words, this class stitch all the code together and contans the necessary setup
    to compute to each individual virus-host pair signals of co-evolution.

    :param virus_directory: Path to the directory of viruses filenames. 1 fasta file = 1 unique virus
    :type virus_directory: str
    :param host_directory: Path to the directory of host filenames. 1 fasta file = 1 unique host
    :type host_directory: str
    :param ext: Extension of fasta filenames. It assumes it is .fasta.
    :type ext: str
    """

    def __init__(
        self, virus_directory, host_directory, ext="fasta", pairs_of_interest=None
    ) -> None:
        self.virus_directory = virus_directory
        self.host_directory = host_directory
        self.ext = ext

        if pairs_of_interest:
            self.pairs_of_interest = pairs_of_interest
        else:
            self.pairs_of_interest = None

        self.features_df = None

    def add_blastn_files(self, blastn_path, spacer_path):
        """Add blastn path for viruses against hosts and viruses against spacers.
        Needed for homology features.

        :param blastn_path: Path to location of blastn output for viruses against hosts
        :type blastn_path: str
        :param spacer_path: Path to location of blastn output for viruses against spacers
        :type blastn_spath: str
        """
        self.blastn_path = blastn_path
        self.spacer_path = spacer_path

    def do_setup(self):
        """This method call on other methods to determine all possible virus-host pairs, get fasta headers,
        read blastn outputs, and computes GC content and k-mer profiles.
        The methods it calls on are defined below.
        """
        print("SETUP - ...indexing fasta filenames for viruses and hosts...")
        self.list_files()

        print("SETUP - ...initialize all pairs...")
        print(self.pairs_of_interest)
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

    def list_files(self):
        """List all fasta file in the virus and host directories."""
        self.virus_filenames = [
            f for f in os.listdir(self.virus_directory) if f.endswith("." + self.ext)
        ]
        self.host_filenames = [
            f for f in os.listdir(self.host_directory) if f.endswith("." + self.ext)
        ]
        self.all_files = self.virus_filenames + self.host_filenames

    def determine_pairs(self):
        """This method determine all possible virus-host pairs. This assume that each virus should be tested against
        each host. The method also initialize the pair dataclass for each possible virus-host pair combination.
        """
        total_interactions = len(self.virus_filenames) * len(self.host_filenames)
        print(f"-------> There are {len(self.virus_filenames)} viral sequences")
        print(f"-------> There are {len(self.host_filenames)} host sequences")
        print(f"-------> Total number of interactions: {total_interactions}")

        # determine all virus-host pair possible (every host is going to be considered for every virus of interest)
        virus_inter = list(
            itertools.chain.from_iterable(
                itertools.repeat(x, len(self.host_filenames))
                for x in self.virus_filenames
            )
        )
        host_inter = self.host_filenames * len(self.virus_filenames)
        pairs = list(zip(virus_inter, host_inter))

        # create list of Pairs
        self.pairs = []
        for pair in list(zip(virus_inter, host_inter)):
            virus = pair[0]
            host = pair[1]
            self.pairs.append(Pairs(virus, host))

    def determine_custom_pairs(self, custom_pairs):
        """The input pair file needs to be virus first then host. Must be separated by tabs."""
        self.pairs = []
        print("reading pairs file")

        with open(custom_pairs, "r") as f:
            lines = [line.rstrip() for line in f]
            for pair in lines:
                split = pair.split("\t")
                virus = split[0]
                host = split[1]

                if virus in self.virus_filenames and host in self.host_filenames:
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
        """Get headers in fasta files. This is because some fasta files might have multiple > in them. This is needed
        to map which blast hit came from which virus/host.
        """
        header_filename = {}

        for virus in self.virus_filenames:
            path = self.virus_directory + virus
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = virus

        for host in self.host_filenames:
            path = self.host_directory + host
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = host

        self.headers = header_filename

    def process_blastn(self):
        """Read blastn file and store hits as a dictionary. Keys of the dictionary are viruses. For each key, then
        there is a list of hosts that had a homology hit with that virus.
        """
        self.blastn = {}

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
        """Read spacer file and store hits as a dictionary. Key of the dictionary are viruses. For each keym then there
        is a list of hosts that had a homology hit with that virus.
        """
        self.spacers = {}

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
                    filter(lambda x: x.startswith(host_partial), self.host_filenames)
                )
                if virus not in self.spacers:
                    self.spacers[virus] = host
                elif host not in self.spacers[virus]:
                    self.spacers[virus].append(host[0])

    def generate_kmer_profiles(self):
        """Method to compute the GC content, kmer profile with k=3, and kmer profile with k=6, for each fasta file
        (viruses and hosts both included).
        """
        self.GCcontent = dict.fromkeys(self.all_files)
        self.k3profiles = dict.fromkeys(self.all_files)
        self.k6profiles = dict.fromkeys(self.all_files)

        for virus in self.virus_filenames:
            path = self.virus_directory + virus
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

        for host in self.host_filenames:
            path = self.host_directory + host
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

    def run_parallel(self, num_procs=6):
        """Call the compute_feature method to compute features for each virus-host pair. This step is parallelized
        to significantly improve run time.

        :param num_procs: Number of core to be used.
        :type num_procs: int
        """
        with Pool(num_procs) as pool:
            results = pool.map(self.compute_feature, self.pairs)

        self.computed_pairs = results

    def compute_feature(self, pair):
        """Calls on feature classes to compute genome level features. Returns the signals of coevolution for
        the given pair.

        :param pair: pair Dataclass.
        :type pair: Pairs
        """
        print(f"-------> current pair: {pair.virus} | {pair.host}")
        pair.homology_hit = self.homology_matches.match(pair.virus, pair.host)

        k3distance = d2Distance(self.k3profiles[pair.virus], self.k3profiles[pair.host])
        k3distance.distance()
        pair.k3dist = k3distance.dist

        k6distance = d2Distance(self.k6profiles[pair.virus], self.k6profiles[pair.host])
        k6distance.distance()
        pair.k6dist = k6distance.dist

        pair.GCdifference = self.GCcontent[pair.virus] - self.GCcontent[pair.host]

        return pair

    def convert_to_dataframe(self):
        """Convert the list of Pairs into a Pandas dataframe."""
        pairs = []

        k3dist = []
        k6dist = []
        GCdiff = []
        Homology = []

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
        self.features_df = self.features_df.set_index("pairs")

    def save_features(self, filename):
        """This is an optional method. This allow user to save the computed features as a tsv file."""
        if self.features_df is None:
            self.convert_to_dataframe()

        self.features_df.to_csv(filename, sep="\t")

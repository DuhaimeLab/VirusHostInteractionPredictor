"""Contains classes to compute genome-wide features.

This module provides:
- KmerProfile: generate k-mer profiles for a give sequence and k-length
- d2Distance: computes the d2* distance between two KmerProfile
- HomologyMatch: check for DNA sequence match between a virus and host from blastn files
"""

from typing import List

import numpy as np


class KmerProfile:
    """Computes the k-mer profile for a given sequence.

    Args:
        seq (str): nucleotide sequence
        k (int): k-length to be used when generating k-mer profile of sequence
    """

    def __init__(self, seq: List[str], k: int) -> None:
        """Initialize class variables."""
        if not seq:
            raise ValueError("seq cannot be an empty list")
        self.k = k
        self.seq = seq
        self.seqlen = sum(len(s) for s in seq)
        self.nucleotides = ["A", "T", "C", "G"]

    def generate_profile(self):
        """Generate k-mer profile for sequence."""
        words = self.generate_kmer_words(self.k)
        self.kmer_words = words
        profile = dict.fromkeys(words, 0)

        for seq in self.seq:
            for i in range(0, len(seq) - self.k + 1):
                kmer = seq[i : i + self.k]
                if kmer in profile:
                    profile[kmer] += 1

        self.profile_counts = np.fromiter(profile.values(), dtype=int)

        # calculate GC content if k = 1 since it"s nucleotide counting
        if self.k == 1:
            AT = self.profile_counts[0] + self.profile_counts[1]
            GC = self.profile_counts[2] + self.profile_counts[3]
            self.GCcontent = GC / (AT + GC) * 100

    def generate_kmer_words(self, length: int) -> list[str]:
        """Generate all possible k-mer words based on parameter length.

        Args:
            length (int): k-length to be used

        Returns:
            list: list of all possible k-mer words
        """
        if length == 1:
            return self.nucleotides
        else:
            result: List[str] = []
            for seq in self.generate_kmer_words(length - 1):
                for base in self.nucleotides:
                    result.append(seq + base)
            return result


class d2Distance:
    """Compute d2* distance between two k-mer profiles.

    The d2* a measurement of how similar or dissimilar those two k-mer profiles are to
    one another, and the value returned ranged from 0 to 1.

    Args:
        seq1_profile (KmerProfile): k-mer profile of sequence 1
        seq2_profile (KmerProfile): k-mer profile of sequence 2
    """

    def __init__(self, seq1_profile: KmerProfile, seq2_profile: KmerProfile) -> None:
        """Initialize class variables."""
        self.k = seq1_profile.k
        self.kmer_words = seq1_profile.kmer_words

        self.seq1_profile = seq1_profile
        self.seq2_profile = seq2_profile

        self.seq1_seqlen = seq1_profile.seqlen
        self.seq2_seqlen = seq2_profile.seqlen

    def distance(self):
        """Compute the distance for the two input k-mer profile.

        Returns:
            float: distance value that is between 0 and 1.
        """
        if self.seq1_profile.k == self.seq2_profile.k:
            self.nucleotide_count()
            self.null()
            self.dist = self.d2star()
            return self.dist

        else:
            print("k-length used to generate the profiles is not consistent")
            self.dist = None
            return self.dist

    def nucleotide_count(self) -> None:
        """Calculate the nucleotide count for sequence 1 and 2. It assumes A, T, C, and G letters."""
        seq1_nuc_count = dict.fromkeys(self.seq1_profile.nucleotides, 0)
        seq2_nuc_count = dict.fromkeys(self.seq2_profile.nucleotides, 0)

        for seq in self.seq1_profile.seq:
            for nuc in seq:
                if nuc in self.seq1_profile.nucleotides:
                    seq1_nuc_count[nuc] += 1
        for seq in self.seq2_profile.seq:
            for nuc in seq:
                if nuc in self.seq2_profile.nucleotides:
                    seq2_nuc_count[nuc] += 1

        # calculate frequencies
        seq1_total = sum(seq1_nuc_count.values(), 0)
        seq2_total = sum(seq2_nuc_count.values(), 0)
        self.seq1_nuc_freq = {
            key: value / seq1_total for key, value in seq1_nuc_count.items()
        }
        self.seq2_nuc_freq = {
            key: value / seq2_total for key, value in seq2_nuc_count.items()
        }

    def null(self) -> None:
        """Computer null expectation of k-mer words based on A, T, C, and G counts."""
        self.nucleotide_count()

        seq1_null = dict.fromkeys(self.kmer_words, 0.0)
        seq2_null = dict.fromkeys(self.kmer_words, 0.0)

        for word in self.kmer_words:
            countA = word.count("A")
            countT = word.count("T")
            countC = word.count("C")
            countG = word.count("G")

            if countA == 0:
                # set to 1 because when computing null, everything is multiplied
                # so 1 would not change the values
                countA = 1
                seq1_nfA = 1
                seq2_nfA = 1
            else:
                seq1_nfA = self.seq1_nuc_freq["A"]
                seq2_nfA = self.seq2_nuc_freq["A"]

            if countT == 0:
                countT = 1
                seq1_nfT = 1
                seq2_nfT = 1
            else:
                seq1_nfT = self.seq1_nuc_freq["T"]
                seq2_nfT = self.seq2_nuc_freq["T"]

            if countC == 0:
                countC = 1
                seq1_nfC = 1
                seq2_nfC = 1
            else:
                seq1_nfC = self.seq1_nuc_freq["C"]
                seq2_nfC = self.seq2_nuc_freq["C"]

            if countG == 0:
                countG = 1
                seq1_nfG = 1
                seq2_nfG = 1
            else:
                seq1_nfG = self.seq1_nuc_freq["G"]
                seq2_nfG = self.seq2_nuc_freq["G"]

            # compute null
            seq1_null[word] = (
                (seq1_nfA**countA)
                * (seq1_nfT**countT)
                * (seq1_nfG**countG)
                * (seq1_nfC**countC)
            )
            seq2_null[word] = (
                (seq2_nfA**countA)
                * (seq2_nfT**countT)
                * (seq2_nfG**countG)
                * (seq2_nfC**countC)
            )

        self.seq1_null_prob = np.fromiter(seq1_null.values(), dtype=float)
        self.seq2_null_prob = np.fromiter(seq2_null.values(), dtype=float)

    def d2star(self) -> float:
        """Normalize the distance.

        Returns:
            float: the distance value normalized, between 0 and 1
        """
        D2star_value = self.D2star()
        numerator = np.sqrt(sum(self.x**2 / self.x_expected)) * np.sqrt(
            sum(self.y**2 / self.y_expected)
        )
        return 0.5 * (1 - (D2star_value / numerator))

    def D2star(self) -> float:
        """Computer distance between two k-mer profiles.

        Returns:
            float: distance measurement
        """
        # calculate expected kmer counts based on null probability multiplied by length of sequence
        self.x_expected = self.seq1_null_prob * self.seq1_profile.seqlen
        self.y_expected = self.seq2_null_prob * self.seq2_profile.seqlen
        # difference in observed kmer counts and expected kmer counts
        self.x = self.seq1_profile.profile_counts - self.x_expected
        self.y = self.seq2_profile.profile_counts - self.y_expected

        # D2* distance between x, the virus, and y, the host
        numerator = self.x * self.y
        denominator = np.sqrt(
            self.seq1_profile.seqlen
            * self.seq2_profile.seqlen
            * self.seq1_null_prob
            * self.seq2_null_prob
        )
        return sum(numerator / denominator)


class HomologyMatch:
    """Class to determine if there are matches in DNA sequences between viruses and hosts of interest using blastn output as input.

    Args:
        virus_host_blastn (dict): dictionary containing viruses as keys and list of host fasta file headers as values.
        virus_spacer_blastn (dict): dictionary containing viruses as keys and list of spacers headers as values.
    """

    def __init__(
        self,
        virus_host_blastn: dict[str, list[str]],
        virus_spacer_blastn: dict[str, list[str]],
    ):
        """Initialize class variables."""
        self.virus_host = virus_host_blastn
        self.virus_spacers = virus_spacer_blastn

    def match(self, virus: str, host: str):
        """Determine if there is a match in homology for given virus-host pair.

        This will check the virus against host blastn results, and the virus against host spacers blastn results.

        Args:
            virus (str): virus filename
            host (str): host filename
        """
        self.check_blastn(virus, host)
        self.check_spacers(virus, host)

        if self.blast or self.spacer:
            self.hit = True
        else:
            self.hit = False

        return self.hit

    def check_blastn(self, virus: str, host: str):
        """Check the virus against host blastn dictionary.

        If there is a match from this blastn dictionary, set self.blast to True.

        Args:
            virus (str): virus filename
            host (str): host filename
        """
        if virus not in self.virus_host:
            self.blast = False
        elif host in self.virus_host[virus]:
            self.blast = True
        else:
            self.blast = False

    def check_spacers(self, virus: str, host: str):
        """Check the virus against spacer blastn dictionary.

        If there is a match from this blastn dictionary, set self.spacer to True.

        Args:
            virus (str): virus filename
            host (str): host filename
        """
        if virus not in self.virus_spacers:
            self.spacer = False
        elif host in self.virus_spacers[virus]:
            self.spacer = True
        else:
            self.spacer = False

'''

'''

import numpy as np


class KmerProfile:

    def __init__(self, seq, k) -> None:
        self.k = k
        self.seq = seq
        self.seqlen = len(seq)
        self.nucleotides = ['A', 'T', 'C', 'G']

        self.profile_counts = None
    
    def generate_profile(self):
        words = self.generate_kmer_words(self.k)
        self.kmer_words = words
        profile = dict.fromkeys(words, 0)

        for i in range (0, len(self.seq) - self.k + 1):
            kmer = self.seq[i: i+self.k]
            if kmer in profile:
                profile[kmer] += 1
        
        self.profile_counts = np.fromiter(profile.values(), dtype=int)
        # code below is not needed
        #total = sum(profile.values(), 0)
        #freqs = {key: value / total for key, value in profile.items()}
        #self.profile_freqs = np.fromiter(freqs.values(), dtype=float)


    def generate_kmer_words(self, length):
        if length == 1:
            return self.nucleotides
        else:
            result = []
            for seq in self.generate_kmer_words(length - 1):
                for base in self.nucleotides:
                    result.append(seq + base)
            return result
        



'''
seq = 'ATCCTAGAGTTAGCCGTAAAAAAAAAAAAAAAA'
test = KmerProfile(seq, 3)
test.generate_profile()
print(test.profile_counts.shape)
print(test.profile_freqs.shape)
'''


class d2Distance:

    def __init__(self, seq1_profile, seq2_profile) -> None:

        # TODO: need check that k-length used in virus profile and host profile match


        self.k = seq1_profile.k
        self.kmer_words = seq1_profile.kmer_words

        self.seq1_profile = seq1_profile
        self.seq2_profile = seq2_profile

        self.seq1_seqlen = seq1_profile.seqlen
        self.seq2_seqlen = seq2_profile.seqlen


    def distance(self):
        self.nucleotide_count()
        self.null()
        self.d2star()
        return self.dist

    
    def nucleotide_count(self) -> None:
        seq1_nuc_count = dict.fromkeys(self.seq1_profile.nucleotides, 0)
        seq2_nuc_count = dict.fromkeys(self.seq2_profile.nucleotides, 0)

        for nuc in self.seq1_profile.seq:
            seq1_nuc_count[nuc] += 1
        for nuc in self.seq2_profile.seq:
            seq2_nuc_count[nuc] += 1

        # calculate frequencies
        seq1_total = sum(seq1_nuc_count.values(), 0)
        seq2_total = sum(seq2_nuc_count.values(), 0)
        self.seq1_nuc_freq = {key: value / seq1_total for key, value in seq1_nuc_count.items()}
        self.seq2_nuc_freq = {key: value / seq2_total for key, value in seq2_nuc_count.items()}


    def null(self):

        self.nucleotide_count()

        seq1_null = dict.fromkeys(self.kmer_words, 0)
        seq2_null = dict.fromkeys(self.kmer_words, 0)

        for word in self.kmer_words:
            
            countA = word.count('A')
            countT = word.count('T')
            countC = word.count('C')
            countG = word.count('G')

            if countA == 0:
                countA = 1
                seq1_nfA = 1
                seq2_nfA = 1
            else:
                seq1_nfA = self.seq1_nuc_freq['A']
                seq2_nfA = self.seq2_nuc_freq['A']
            
            if countT == 0:
                countT = 1
                seq1_nfT = 1
                seq2_nfT = 1
            else:
                seq1_nfT = self.seq1_nuc_freq['T']
                seq2_nfT = self.seq2_nuc_freq['T']

            if countC == 0:
                countC = 1
                seq1_nfC = 1
                seq2_nfC = 1
            else:
                seq1_nfC = self.seq1_nuc_freq['C']
                seq2_nfC = self.seq2_nuc_freq['C']

            if countG == 0:
                countG = 1
                seq1_nfG = 1
                seq2_nfG = 1
            else:
                seq1_nfG = self.seq1_nuc_freq['G']
                seq2_nfG = self.seq2_nuc_freq['G']
            
            # compute null 
            seq1_null[word] = (seq1_nfA ** countA) * (seq1_nfT ** countT) * (seq1_nfG ** countG) * (seq1_nfC ** countC)
            seq2_null[word] = (seq2_nfA ** countA) * (seq2_nfT ** countT) * (seq2_nfG ** countG) * (seq2_nfC ** countC)

        self.seq1_null_prob = np.fromiter(seq1_null.values(), dtype=float)
        self.seq2_null_prob = np.fromiter(seq2_null.values(), dtype=float)


    def d2star(self):

        D2star_value = self.D2star()
        numerator = np.sqrt(sum(self.x ** 2 / self.x_expected)) * np.sqrt(sum(self.y ** 2 / self.y_expected))
        self.dist = 0.5 * (1 - (D2star_value / numerator))


    def D2star(self):

        # difference in observed kmer counts and expected kmer counts
        self.x_expected = self.seq1_null_prob * self.seq1_profile.seqlen
        self.y_expected = self.seq2_null_prob * self.seq2_profile.seqlen

        self.x = self.seq1_profile.profile_counts - self.x_expected
        self.y = self.seq2_profile.profile_counts - self.y_expected

        # D2* distance between x, the virus, and y, the host
        numerator = self.x * self.y
        denominator = np.sqrt(self.seq1_profile.seqlen * self.seq2_profile.seqlen * self.seq1_null_prob * self.seq2_null_prob)
        return sum(numerator / denominator)






'''
virus_seq = 'GGGCCCCCTTTAAAA'
host_seq = 'AAATTTCCCGGG'

virus_kmer_profile = KmerProfile(virus_seq, 6)
virus_kmer_profile.generate_profile()
host_kmer_profile = KmerProfile(host_seq, 6)
host_kmer_profile.generate_profile()

test = d2Distance(virus_kmer_profile, host_kmer_profile)
test.distance()
print(test.distance())

test2 = d2Distance(host_kmer_profile, virus_kmer_profile)
test2.distance
print(test2.distance())


virus_kmer_profile = KmerProfile(virus_seq, 1)
virus_kmer_profile.generate_profile()
print(virus_kmer_profile.nucleotides)
print(virus_kmer_profile.profile_counts)

'''


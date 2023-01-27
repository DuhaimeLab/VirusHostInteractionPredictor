'''

'''

import numpy as np


class KmerProfile:

    def __init__(self, seq, k) -> None:
        self.k = k
        self.seq = seq
        self.seqlen = len(seq)
        self.nucleotides = ['A', 'T', 'C', 'G']

    
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

    def __init__(self, virus_profile, host_profile) -> None:

        # TODO: need check that k-length used in virus profile and host profile match


        self.k = virus_profile.k
        self.kmer_words = virus_profile.kmer_words

        self.virus_profile = virus_profile
        self.host_profile = host_profile

        self.virus_seqlen = virus_profile.seqlen
        self.host_seqlen = host_profile.seqlen


    def distance(self):
        self.nucleotide_count()
        self.null()
        return self.d2star()

    
    def nucleotide_count(self) -> None:
        virus_nuc_count = dict.fromkeys(self.virus_profile.nucleotides, 0)
        host_nuc_count = dict.fromkeys(self.host_profile.nucleotides, 0)

        for nuc in self.virus_profile.seq:
            virus_nuc_count[nuc] += 1
        for nuc in self.host_profile.seq:
            host_nuc_count[nuc] += 1

        # calculate frequencies
        virus_total = sum(virus_nuc_count.values(), 0)
        host_total = sum(host_nuc_count.values(), 0)
        self.virus_nuc_freq = {key: value / virus_total for key, value in virus_nuc_count.items()}
        self.host_nuc_freq = {key: value / host_total for key, value in host_nuc_count.items()}


    def null(self):

        self.nucleotide_count()

        virus_null = dict.fromkeys(self.kmer_words, 0)
        host_null = dict.fromkeys(self.kmer_words, 0)

        for word in self.kmer_words:
            
            countA = word.count('A')
            countT = word.count('T')
            countC = word.count('C')
            countG = word.count('G')

            if countA == 0:
                countA = 1
                virus_nfA = 1
                host_nfA = 1
            else:
                virus_nfA = self.virus_nuc_freq['A']
                host_nfA = self.host_nuc_freq['A']
            
            if countT == 0:
                countT = 1
                virus_nfT = 1
                host_nfT = 1
            else:
                virus_nfT = self.virus_nuc_freq['T']
                host_nfT = self.host_nuc_freq['T']

            if countC == 0:
                countC = 1
                virus_nfC = 1
                host_nfC = 1
            else:
                virus_nfC = self.virus_nuc_freq['C']
                host_nfC = self.host_nuc_freq['C']

            if countG == 0:
                countG = 1
                virus_nfG = 1
                host_nfG = 1
            else:
                virus_nfG = self.virus_nuc_freq['G']
                host_nfG = self.host_nuc_freq['G']
            
            # compute null 
            virus_null[word] = (virus_nfA ** countA) * (virus_nfT ** countT) * (virus_nfG ** countG) * (virus_nfC ** countC)
            host_null[word] = (host_nfA ** countA) * (host_nfT ** countT) * (host_nfG ** countG) * (host_nfC ** countC)

        self.virus_null_prob = np.fromiter(virus_null.values(), dtype=float)
        self.host_null_prob = np.fromiter(host_null.values(), dtype=float)


    def d2star(self):

        D2star_value = self.D2star()
        numerator = np.sqrt(sum(self.x ** 2 / self.x_expected)) * np.sqrt(sum(self.y ** 2 / self.y_expected))
        return 0.5 * (1 - (D2star_value / numerator))


    def D2star(self):

        # difference in observed kmer counts and expected kmer counts
        self.x_expected = self.virus_null_prob * self.virus_profile.seqlen
        self.y_expected = self.host_null_prob * self.host_profile.seqlen

        self.x = self.virus_profile.profile_counts - self.x_expected
        self.y = self.host_profile.profile_counts - self.y_expected

        # D2* distance between x, the virus, and y, the host
        numerator = self.x * self.y
        denominator = np.sqrt(self.virus_profile.seqlen * self.host_profile.seqlen * self.virus_null_prob * self.host_null_prob)
        return sum(numerator / denominator)







virus_seq = 'GGGCCCTTTAAA'
host_seq = 'AAATTTCCCGGG'

virus_kmer_profile = KmerProfile(virus_seq, 6)
virus_kmer_profile.generate_profile()
host_kmer_profile = KmerProfile(host_seq, 6)
host_kmer_profile.generate_profile()

test = d2Distance(virus_kmer_profile, host_kmer_profile)
test.distance()
print(test.distance())



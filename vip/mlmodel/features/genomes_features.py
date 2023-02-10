'''

'''

import numpy as np


class KmerProfile:
    '''
    The KmerProfile class takes for input a sequence and a k-mer length. Its purpose
    is to compute the k-mer profile for the given sequence and k length. 
    
    :param seq: A sequence of DNA from which the k-mer profile is to be calculated
    :type seq: str
    :param k: k-length to use to compute the k-mers
    :type k: int
     '''

    def __init__(self, seq, k) -> None:
        self.k = k
        self.seq = seq
        self.seqlen = len(seq)
        self.nucleotides = ['A', 'T', 'C', 'G']

        self.profile_counts = None
    
    def generate_profile(self):
        '''
        Method that generates the k-mer profile for a sequence. No parameters
        are needed since they are stored in self when initializing object. 

        If the k-length is set to 1, then this method also calculate GC content.

        '''

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

        # calculate GC content if k = 1 since it's nucleotide counting
        if self.k == 1:
            AT = self.profile_counts[0] + self.profile_counts[1]
            GC = self.profile_counts[2] + self.profile_counts[3]
            self.GCcontent = GC / (AT + GC) * 100 


    def generate_kmer_words(self, length):
        '''
        Function to generate all possible k-mer words possible for a given length

        :param length: int, k-length to be used to generate the words
        '''
        if length == 1:
            return self.nucleotides
        else:
            result = []
            for seq in self.generate_kmer_words(length - 1):
                for base in self.nucleotides:
                    result.append(seq + base)
            return result
        

class d2Distance:
    '''
    This class allows user to compute the d2* distance between two k-mer profiles of interest. 
    The d2* distance is a measurement of how similar or dissimilar those two k-mer profiles are to 
    one another, and the value returned ranged from 0 to 1. 
    
    :param seq1_profile: k-mer profile of first sequence
    :type seq1_profile: class:KmerProfile
    :param seq2_profile: k-mer profile of second sequence (to compare with first sequence)
    :type seq2_profile: class:KmerProfile
    '''

    def __init__(self, seq1_profile, seq2_profile) -> None:

        # TODO: need check that k-length used in virus profile and host profile match


        self.k = seq1_profile.k
        self.kmer_words = seq1_profile.kmer_words

        self.seq1_profile = seq1_profile
        self.seq2_profile = seq2_profile

        self.seq1_seqlen = seq1_profile.seqlen
        self.seq2_seqlen = seq2_profile.seqlen


    def distance(self):
        '''
        Method to call all the other methods. It first computes the nucleotide_count, which is needed
        to generate the null expectation of k-mer profile. This, in turn, is needed to compute the 
        d2* distance

        :return: float, between 0 and 1.
        '''
        self.nucleotide_count()
        self.null()
        self.d2star()
        return self.dist

    
    def nucleotide_count(self) -> None:
        '''
        Determine the nucleotide count for sequence 1 and 2. It assumes A, T, C, and G letters. 
        '''
        seq1_nuc_count = dict.fromkeys(self.seq1_profile.nucleotides, 0)
        seq2_nuc_count = dict.fromkeys(self.seq2_profile.nucleotides, 0)

        for nuc in self.seq1_profile.seq:
            if nuc in self.seq1_profile.nucleotides: 
                seq1_nuc_count[nuc] += 1
        for nuc in self.seq2_profile.seq:
            if nuc in self.seq2_profile.nucleotides: 
                seq2_nuc_count[nuc] += 1

        # calculate frequencies
        seq1_total = sum(seq1_nuc_count.values(), 0)
        seq2_total = sum(seq2_nuc_count.values(), 0)
        self.seq1_nuc_freq = {key: value / seq1_total for key, value in seq1_nuc_count.items()}
        self.seq2_nuc_freq = {key: value / seq2_total for key, value in seq2_nuc_count.items()}


    def null(self):
        '''
        Compute the null expectation of k-mer words based on A, T, C, and G counts
        '''

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
        '''
        Normalize the distance between two k-mer profile
        '''
        D2star_value = self.D2star()
        numerator = np.sqrt(sum(self.x ** 2 / self.x_expected)) * np.sqrt(sum(self.y ** 2 / self.y_expected))
        self.dist = 0.5 * (1 - (D2star_value / numerator))


    def D2star(self):
        '''
        Compute the distance between two kmer profiles
        '''
        # calculate expected kmer counts based on null probability multiplied by length of sequence
        self.x_expected = self.seq1_null_prob * self.seq1_profile.seqlen
        self.y_expected = self.seq2_null_prob * self.seq2_profile.seqlen
        # difference in observed kmer counts and expected kmer counts
        self.x = self.seq1_profile.profile_counts - self.x_expected
        self.y = self.seq2_profile.profile_counts - self.y_expected

        # D2* distance between x, the virus, and y, the host
        numerator = self.x * self.y
        denominator = np.sqrt(self.seq1_profile.seqlen * self.seq2_profile.seqlen * self.seq1_null_prob * self.seq2_null_prob)
        return sum(numerator / denominator)



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


print(test.x_expected)



class HomologyMatch:

    def __init__(self):
        pass

    def read_blastn():
        pass

    def read_spacers():
        pass

    def check_match():
        pass


test = HomologyMatch()

#test.read_blastn()
#test.read_spacers()

#test.check_match(virus_filename, host_filename)

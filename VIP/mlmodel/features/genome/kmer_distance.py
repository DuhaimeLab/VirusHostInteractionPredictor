'''

'''

def generate_kmers(length):
    if length == 1:
        return ['A', 'T', 'C', 'G']
    else:
        result = []
        for seq in generate_kmers(length - 1):
            for base in ['A', 'T', 'C', 'G']:
                result.append(seq + base)
        return result



def kmer_profile(seq, k, freq=False):

    words = generate_kmers(k)
    profile = dict.fromkeys(words, 0)

    for i in range(0, len(seq) - k + 1): 
        kmer = seq[i: i+k]
        if kmer in profile:
            profile[kmer] += 1
        

# function below is for cases where sequence is divided into smaller separate contigs
def combine_kmer_profile():
    pass

'''

'''


# calculate counts of A, T, C, and G
def nucleotide_counts(seq: str) -> dict:

    # initialize dictionary to keep track of nucleotide counts
    nucleotides = ['A', 'T', 'C', 'G']
    counts = dict.fromkeys(nucleotides, 0)

    for nuc in seq:
        if nuc in nucleotides:
            counts[nuc] += 1
    
    return counts


# compute GC content from fasta file
def GC_content(seq: str) -> float:

    counts = nucleotide_counts(seq)

    Acount = counts['A']
    Tcount = counts['T']
    Gcount = counts['G']
    Ccount = counts['C']

    content = ((Gcount + Ccount) / (Acount + Tcount + Gcount + Ccount)) * 100 

    return content
    

# return difference between virus and host GC content
def GC_difference(GCvirus: float, GChost: float) -> float:
    diff = GCvirus - GChost
    return diff


'''
test_seq = 'ATCCTAGGATCTAGAHDAADTDGTTCGA'
test_seq2 = 'AATC'
print(GC_content(test_seq))
print(GC_content(test_seq2))
'''

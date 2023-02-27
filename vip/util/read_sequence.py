from Bio import SeqIO

import os
print(os.getcwd())

def read_sequence(path):
    '''
    Return the sequence as a string. 
    Take path where file reside as input. 
    '''
    for record in SeqIO.parse(path, 'fasta'):
        return str(record.seq)
    

def read_headers(path):
    '''
    '''
    split = path.split('/')
    filename = split[-1]
    result = [filename]
    for record in SeqIO.parse(path, 'fasta'):
        result.append(record.id)
    return result



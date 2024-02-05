'''Test suite for the ComputeFeatures module.'''

from vhip.mlmodel.compute_ml_features import ComputeFeatures, Pairs
from vhip.mlmodel.genomes_features import KmerProfile

import pandas as pd

test_virus_directory = 'tests/datatests/sequences/virus_seqs/'
test_host_directory = 'tests/datatests/sequences/host_seqs/'


def test_ComputeFeatures_list_files():
    '''Test that listed files are correct.'''
    all_filenames = [
        'GCA_003931015.1_ASM393101v1_genomic.fasta',
        'GCA_003927235.1_ASM392723v1_genomic.fasta',
        'GCA_003931415.1_ASM393141v1_genomic.fasta',
        'GCA_003344205.1_ASM334420v1_genomic.fasta',
        'GCA_005146815.1_ASM514681v1_genomic.fna.fasta',
        'GCA_001974575.1_ASM197457v1_genomic.fna.fasta',
        'GCA_002875995.1_ASM287599v1_genomic.fna.fasta'
        ]

    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    assert len(test.all_files) == 7
    assert len(test.virus_filenames) == 4
    assert len(test.host_filenames) == 3
    assert test.all_files == all_filenames


def test_ComputeFeatures_pairs():
    '''Test all pairs are generated correctly.'''
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.determine_pairs()
    assert len(test.pairs) == 12
    assert all(isinstance(pair, Pairs) for pair in test.pairs)


def test_ComputeFeatures_get_headers():
    '''Test that headers are retrieved correctly.'''
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.get_headers()

    assert test.headers['MG592522.1'] == 'GCA_003344205.1_ASM334420v1_genomic.fasta'
    assert test.headers['MKKP01000001.1'] == 'GCA_001974575.1_ASM197457v1_genomic.fna.fasta'
    assert test.headers['MKKP01000002.1'] == 'GCA_001974575.1_ASM197457v1_genomic.fna.fasta'
    assert test.headers['MKKP01000003.1'] == 'GCA_001974575.1_ASM197457v1_genomic.fna.fasta'


def test_ComputeFeatures_add_blastn_files():
    '''Test that blastn files are added correctly.'''
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    blastn_filename = 'tests/datatests/blastn.tsv'
    spacer_filename = 'tests/datatests/spacer.tsv'
    test.add_blastn_files(blastn_filename, spacer_filename)
    assert test.blastn_path == blastn_filename
    assert test.spacer_path == spacer_filename


def test_ComputeFeatures_process_blastn():
    '''Test that blastn file is processed correctly.'''
    expected_results = {
        'GCA_003344205.1_ASM334420v1_genomic.fasta': ['GCA_003931015.1_ASM393101v1_genomic.fasta'],
        'GCA_003931015.1_ASM393101v1_genomic.fasta': ['GCA_005146815.1_ASM514681v1_genomic.fna.fasta', 'GCA_003931015.1_ASM393101v1_genomic.fasta']}
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.get_headers()
    test.add_blastn_files('tests/datatests/blastn_phagevhost.tsv', '')
    test.process_blastn()
    assert test.blastn == expected_results


def test_ComputeFeatures_process_spacers():
    '''Test that spacers are processed correctly.'''
    expected_results = {
        'GCA_003344205.1_ASM334420v1_genomic.fasta': ['GCA_002875995.1_ASM287599v1_genomic.fna.fasta',
                                                      'GCA_002875995.1_ASM287599v1_genomic.fna.fasta',
                                                      'GCA_002875995.1_ASM287599v1_genomic.fna.fasta'
                                                      ],
        'GCA_003931015.1_ASM393101v1_genomic.fasta': ['GCA_002875995.1_ASM287599v1_genomic.fna.fasta',
                                                      'GCA_002875995.1_ASM287599v1_genomic.fna.fasta'
                                                      ]
    }
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.get_headers()
    test.add_blastn_files('', 'tests/datatests/blastn_phagevspacer.tsv')
    test.process_spacers()
    assert test.spacers == expected_results


def test_ComputeFeatures_generate_kmer_profiles():
    '''Test the kmer profiles are properly generated.'''
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.determine_pairs()
    test.generate_kmer_profiles()
    assert isinstance(test.k6profiles, dict)
    assert isinstance(test.k6profiles['GCA_003344205.1_ASM334420v1_genomic.fasta'], KmerProfile)


def test_ComputeFeatures_complete_pipeline():
    '''Check the complete pipeline for ComputeFeatures is working as intended.'''
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.add_blastn_files('tests/datatests/blastn_phagevhost.tsv', 'tests/datatests/blastn_phagevspacer.tsv')
    test.do_setup()
    test.run_parallel()
    test.convert_to_dataframe()
    assert isinstance(test.features_df, pd.DataFrame)
    assert test.features_df.shape[0] == 12
    assert test.features_df.shape[1] == 4



if __name__ == '__main__':
    test = ComputeFeatures(test_virus_directory, test_host_directory)
    test.list_files()
    test.determine_pairs()
    test.generate_kmer_profiles()
    print(test.k6profiles['GCA_003344205.1_ASM334420v1_genomic.fasta'])

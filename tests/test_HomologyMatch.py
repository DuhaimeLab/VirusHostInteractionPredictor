'''Pytest for HomologyMatch module.'''

from vhip.mlmodel.genomes_features import HomologyMatch

def test_HomologyMatch_check_blastn():
    '''Test if match is retrieved successfully.'''
    virus_host_blastn_test = {'NC_1': ['host1', 'host5'], 'NC_2': ['host1', 'host2']}
    virus_host_spacer_test = {'NC_1': ['spacer10'], 'NC_2': ['spacer1']}

    match_object = HomologyMatch(virus_host_blastn_test, virus_host_spacer_test)
    assert match_object.match('NC_1', 'host1') is True
    assert match_object.match('NC_1', 'host2') is False
    assert match_object.match('NC_1', 'host7') is False
    assert match_object.match('NC_donotexsist', 'host1') is False


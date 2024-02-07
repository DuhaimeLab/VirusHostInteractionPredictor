'''Test suite for the PredictInteractions module.'''

from vhip.predict_interactions import PredictInteractions

MODEL = 'tests/datatests/gbrt.skops'
test_virus_directory = 'tests/datatests/sequences/virus_seqs/'
test_host_directory = 'tests/datatests/sequences/host_seqs/'


def test_PredictInteractions_complete_pipeline():
    '''Test that the complete pipeline to compute and predict virus-host interaction is working.'''
    test = PredictInteractions(test_virus_directory, test_host_directory)
    test.load_model(MODEL)
    test.add_blastn_files('tests/datatests/blastn_phagevhost.tsv', 'tests/datatests/blastn_phagevspacer.tsv')
    test.do_setup()
    test.run_parallel()
    test.predict()

    assert test.predictions is not None
    assert test.scores is not None


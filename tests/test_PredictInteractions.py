"""Tests for the PredictInteractions module."""

from vhip.mlmodel.build import BuildModel
from vhip.predict_interactions import PredictInteractions

test_virus_genome_dir = "tests/datatests/sequences/virus_genomes/"
test_host_genome_dir = "tests/datatests/sequences/host_genomes/"
ml_training = "tests/datatests/ml_input.csv"


def test_PredictInteractions_complete_pipeline():
    """Test that the complete pipeline to compute and predict virus-host interaction is working."""
    model = BuildModel(ml_training)
    VHIP = model.build()

    test = PredictInteractions(test_virus_genome_dir, test_host_genome_dir)
    test.model = VHIP
    test.add_blastn_files(
        "tests/datatests/blastn_phagevhost.tsv",
        "tests/datatests/blastn_phagevspacer.tsv",
    )
    test.do_setup()
    test.run_parallel()
    test.predict()

    assert test.predictions is not None  # pyright: ignore
    assert test.scores is not None  # pyright: ignore

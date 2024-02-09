'''Test suite for the BuildModel module.'''

from vhip.mlmodel.build import BuildModel
from sklearn.ensemble import GradientBoostingClassifier


def test_BuildModel_init():
    '''Test that the BuildModel class is initialized correctly.'''
    test = BuildModel('tests/datatests/ml_input.csv')

    # check the parameters are set correctly
    assert test.default_max_depth == 15
    assert test.default_learning_rate == 0.75
    assert test.default_loss == "exponential"

    # check dataframe is formatted correctly
    assert test.ml_input.shape == (6925, 4) #pyright: ignore
    assert test.ml_target.shape == (6925,) #pyright: ignore


def test_BuildModel_build():
    '''Test that the BuildModel class is built correctly.'''
    test = BuildModel('tests/datatests/ml_input.csv')
    test.build()

    assert isinstance(test.gbrt, GradientBoostingClassifier) #pyright: ignore

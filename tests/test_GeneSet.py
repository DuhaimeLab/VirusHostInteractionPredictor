"""Pytest for GeneSet methods in gene_features module."""

import pytest

from vhip.mlmodel.gene_features import GeneSet


def test_GeneSet_init():
    """Test code to create GeneSet object and initialize class attributes."""
    # test 1 - test GeneSet object creation generates expected number of genes in genes attribute
    test_GeneSet = GeneSet("tests/datatests/test_annotated_genes.ffn")
    assert len(test_GeneSet.genes) == 3

    # test 2 - test GeneSet object creation generates expected gene attributes for a gene
    assert (
        test_GeneSet.genes[2].seq
        == "TTGGTTGAAGAAGTAGTTGTAGATGGCGACATCACATTAGGACAATTTCTAAAGACGGAAGGTATTATCGAATCTGGCGGGCAAGCGAAATGGTTCTTAAATGAGTTTGAAGTATTGTTAAACAATACGCGTGAAACACGCCGTGGTAAAAAGTTAAGCCATCGTGACACAATTGAGATACCAGAAATACCTGAAGTGGGTTCATTTGTGATTTTGCATCAAGGTGAAGAATGA"
    )
    assert test_GeneSet.genes[2].gene_id == "ABDEAL_00015"
    assert test_GeneSet.genes[2].gene_product == "S4 domain-containing protein YaaA"

    # test 3 - test Exception is raised if gene file is empty
    with pytest.raises(Exception):
        GeneSet("tests/datatests/test_empty_file.ffn")

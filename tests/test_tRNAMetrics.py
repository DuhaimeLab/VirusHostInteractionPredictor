"""Pytest for tRNAMetrics methods in gene_features module."""

import math

import scipy  # pyright: ignore[reportMissingTypeStubs]

from vhip.mlmodel.gene_features import GeneSet, tRNAMetrics

virus_geneset = GeneSet("tests/datatests/test_virus_tRNA_genes.ffn")
host_geneset = GeneSet("tests/datatests/test_host_tRNA_genes.ffn")

def test_tRNAMetrics_init():
    """Test code to create tRNAMetrics object and initialize class attributes."""
    test_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    assert hasattr(test_tRNAMetrics, "virus_GeneSet")
    assert hasattr(test_tRNAMetrics, "host_GeneSet")
    assert hasattr(test_tRNAMetrics.virus_GeneSet, "tRNA_frq_aa")
    assert hasattr(test_tRNAMetrics.host_GeneSet, "tRNA_frq_aa")

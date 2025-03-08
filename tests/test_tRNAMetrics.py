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

def test_tRNAMetrics_virus_TAAI_defaults():
    """Test code to calculate virus amino acid accordance with tRNA availability, using default setting (include virus tRNA genes)."""
    test_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test_tRNAMetrics.virus_TAAI()

    # test 1: check that amino acid frequency has been calculated for virus GeneSet attribute
    assert hasattr(test_tRNAMetrics.virus_GeneSet, "aa_frq")

    # test 2: check TAAI logic (spearman rank between identical data should be 1, and between opposite data should be -1)
    data = test_tRNAMetrics.host_GeneSet.tRNA_frq_aa.values()
    res1 = scipy.stats.spearmanr(list(data), list(data))
    assert math.isclose(res1.statistic, 1.0, rel_tol=1e-6)
    res2 = scipy.stats.spearmanr(list(data), list(data)[::-1])
    assert math.isclose(res2.statistic, -1.0, rel_tol=1e-6)

    # test 3: check correct correlation coefficient between virus amino acid frequency and host tRNA frequency
    assert math.isclose(test_tRNAMetrics.virusTAAI_hosttRNA, 0.7814187052403264, rel_tol=1e-6)

    # test 4: check correc correlation coefficient between virus amino acid frequency and TOTAL tRNA frequency
    assert math.isclose(test_tRNAMetrics.virusTAAI_totaltRNA, 0.7814187052403264, rel_tol=1e-6)



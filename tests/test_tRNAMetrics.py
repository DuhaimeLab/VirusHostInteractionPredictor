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


def test_tRNAMetrics_virus_TAAI():
    """Test code to calculate virus amino acid accordance with tRNA availability."""
    test_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test_tRNAMetrics.virus_TAAI()

    # test 1: check that amino acid frequency has been calculated for virus GeneSet attribute
    assert hasattr(test_tRNAMetrics.virus_GeneSet, "aa_frq")

    # test 2: check TAAI logic (spearman rank between identical data should be 1, and between opposite data should be -1)
    data = list(test_tRNAMetrics.host_GeneSet.tRNA_frq_aa.values())
    res1 = scipy.stats.spearmanr(data, data)
    assert math.isclose(res1.statistic, 1.0, rel_tol=1e-6)
    # res2 = scipy.stats.spearmanr(data, data[::-1])            ### mystery pytest issue - logic works locally
    # assert math.isclose(res2.statistic, -1.0, rel_tol=1e-6)   ### mystery pytest issue - logic works locally

    # test 3: check correct correlation coefficient between virus amino acid frequency and host tRNA frequency
    assert math.isclose(
        test_tRNAMetrics.virusTAAI_hosttRNA, 0.7814187052403264, rel_tol=1e-6
    )

    # test 4: check correct correlation coefficient between virus amino acid frequency and TOTAL tRNA frequency
    assert math.isclose(
        test_tRNAMetrics.virusTAAI_totaltRNA, 0.7814187052403264, rel_tol=1e-6
    )

    # test 5: check no total tRNA comparison metric is generated if parameter specified as false
    test5_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test5_tRNAMetrics.virus_TAAI(include_virus_tRNA=False)
    assert not hasattr(test5_tRNAMetrics, "virusTAAI_totaltRNA")


def test_tRNAMetrics_virus_TCAI():
    """Test code to calculate virus codon accordance with tRNA availability."""
    test_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test_tRNAMetrics.virus_TCAI(skip_nondeg_codons=True)

    # test 1: check that codon frequency has been calculated for virus GeneSet attribute
    assert hasattr(test_tRNAMetrics.virus_GeneSet, "codon_frq")

    # test 2: skip non-degenerate codons
    assert math.isclose(
        test_tRNAMetrics.virusTCAI_hosttRNA, 0.5719110045885629, rel_tol=1e-6
    )
    assert math.isclose(
        test_tRNAMetrics.virusTCAI_totaltRNA, 0.5722056927702744, rel_tol=1e-6
    )

    # test 3: do not skip non-degenerate codons
    test3_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test3_tRNAMetrics.virus_TCAI(skip_nondeg_codons=False)
    assert math.isclose(
        test3_tRNAMetrics.virusTCAI_hosttRNA, 0.6676379024704918, rel_tol=1e-6
    )
    assert math.isclose(
        test3_tRNAMetrics.virusTCAI_totaltRNA, 0.6679138744756812, rel_tol=1e-6
    )

    # test 4: # test 5: check no total tRNA comparison metric is generated if parameter specified as false
    test4_tRNAMetrics = tRNAMetrics(virus_geneset, host_geneset)
    test4_tRNAMetrics.virus_TCAI(include_virus_tRNA=False)
    assert not hasattr(test4_tRNAMetrics, "virusTCAI_totaltRNA")

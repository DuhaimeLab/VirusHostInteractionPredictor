"""Pytest for tRNAMetrics methods in gene_features module."""

import math

import scipy  # pyright: ignore[reportMissingTypeStubs]

from vhip.mlmodel.gene_features import GeneSet, tRNAMetrics

# Set up test data
host_GS = GeneSet("tests/datatests/test_tRNA_genes.json")
virus_GS = GeneSet("tests/datatests/test_short_genes.json")

virus_GS.amino_acid_frequency()
virus_GS.tRNA_counts()
host_GS.tRNA_counts()

virus_aa_frq = virus_GS.aa_frq
host_tRNA_dict_aa = host_GS.tRNA_dict_aa_1_letter
virus_tRNA_dict_aa = virus_GS.tRNA_dict_aa_1_letter

def test_tRNAMetrics_init():
    """Test code to create tRNAMetrics object and initialize class attributes."""

def test_tRNAMetrics_virus_TAAI():
    """Test code to calculate virus amino acid accordance with tRNA availability."""
    test_tRNAMetrics = tRNAMetrics()
    test_tRNAMetrics.virus_TAAI(virus_aa_frq=virus_aa_frq, host_tRNA_dict_aa=host_tRNA_dict_aa, virus_tRNA_dict_aa=virus_tRNA_dict_aa)

    # test 1: check TAAI logic (spearman rank between identical data should be 1, and between opposite data should be -1)
    sorted_keys = sorted(virus_aa_frq)
    data = [host_tRNA_dict_aa[key] for key in sorted_keys]
    res1 = scipy.stats.spearmanr(data, data)
    assert math.isclose(res1.statistic, 1.0, rel_tol=1e-6)

    # test 2: check correct correlation coefficient between virus amino acid frequency and host tRNA frequency
    assert math.isclose(
        test_tRNAMetrics.virusTAAI_hosttRNA, -0.1972026594366539, rel_tol=1e-6
    )

    # test 3: check correct correlation coefficient between virus amino acid frequency and TOTAL tRNA frequency
    assert math.isclose(
        test_tRNAMetrics.virusTAAI_totaltRNA, 0.07562376774775251, rel_tol=1e-6
    )

    # test 4: check no total tRNA comparison metric is generated if no virus tRNA dict provided
    test4_tRNAMetrics = tRNAMetrics()
    test4_tRNAMetrics.virus_TAAI(virus_aa_frq=virus_aa_frq, host_tRNA_dict_aa=host_tRNA_dict_aa)
    assert not hasattr(test4_tRNAMetrics, "virusTAAI_totaltRNA")











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

import pytest

from goreverselookup.pvalcalc import pvalue_calculate


def test_available_pvalue_calculate():
    assert (
        pytest.approx(pvalue_calculate(1, 2, 3, 40, "fisher_scipy_stats"))
        == 0.1461538461538462
    )
    assert (
        pytest.approx(pvalue_calculate(1, 2, 3, 40, "binomial_scipy_stats"))
        == 0.1443750
    )


def test_pvalue_calculate_exception():
    with pytest.raises(ValueError):
        pvalue_calculate(1, 1, 1, 1, "notamethod")


def test_greater_alternative_fisher_scipy():
    pvals = []
    for i in [0, 1, 2, 5, 10, 20]:
        pvalue = pvalue_calculate(i, 20, 100, 1000, "fisher_scipy_stats")
        print(pvalue)
        pvals.append(pvalue)
    assert all(pvals[i] < pvals[i - 1] for i in range(1, len(pvals)))

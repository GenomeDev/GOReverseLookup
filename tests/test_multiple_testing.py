import pytest

from goreverselookup.multiple_testing import multiple_correction


@pytest.fixture
def test_pvals():
    return [0.1, 0.2, 0.5, 0.05, 0.01]


@pytest.mark.parametrize(
    "method",
    [
        "bonferroni",
        "sm_bonferroni",
        "sm_sidak",
        "sm_holm-sidak",
        "sm_holm",
        "sm_simes-hochberg",
        "sm_hommel",
        "sm_fdr_bh",
        "sm_fdr_by",
        "sm_fdr_tsbh",
        "sm_fdr_tsbky",
    ],
)
def test_available_multiple_correction(test_pvals, method):
    corrected_pvals = multiple_correction(test_pvals, method)
    assert len(corrected_pvals) == 5


def test_exceptions_multiple_correction():
    for method in ["sm_notthere", "notin"]:
        with pytest.raises(ValueError):
            multiple_correction([1], method)
    for method in ["holm", "sidak", "fdr_bh"]:
        with pytest.raises(NotImplementedError):
            multiple_correction([1], method)

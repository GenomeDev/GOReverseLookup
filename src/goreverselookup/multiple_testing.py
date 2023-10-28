from statsmodels.stats.multitest import (  # TODO: in future only import when needed and only once
    multipletests,
)


def bonferroni(pvals: list[float], a) -> list[float]:
    """bonferroni correction

    Args:
        pvals (_type_): _description_
        a (_type_): _description_

    Returns:
        _type_: _description_
    """
    n = len(pvals)
    corrected_pvals = [min(p * n, 1) for p in pvals]

    return corrected_pvals


def multiple_correction(pvals: list[float], method: str, a=0.05) -> list[float]:
    """_summary_

    Args:
        pvals (_type_): _description_
        method (str): selected method. statsmodels method are prefixed by sm_
        a (float, optional): _description_. Defaults to 0.05.

    Raises:
        NotImplemented: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    corrected_pvals: list[float]
    if method == "bonferroni":
        corrected_pvals = bonferroni(pvals, a)
    elif method == "holm":
        raise NotImplementedError
    elif method == "sidak":
        raise NotImplementedError
    elif method == "fdr_bh":
        raise NotImplementedError
    elif "sm_" in method:
        test_name = method[3:]
        if test_name not in [
            "bonferroni",
            "sidak",
            "holm-sidak",
            "holm",
            "simes-hochberg",
            "hommel",
            "fdr_bh",
            "fdr_by",
            "fdr_tsbh",
            "fdr_tsbky",
        ]:
            raise ValueError(f"{method} not in statsmodels multipletests")
        corrected_pvals = list(multipletests(pvals, a, "bonferroni")[1])
    else:
        raise ValueError(f"{method} not in available methods")

    return corrected_pvals

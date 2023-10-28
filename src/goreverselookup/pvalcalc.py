"""Options for calculating uncorrected p-values."""
from scipy import stats


def fisherscipystats(study_count, study_n, pop_count, pop_n) -> float:
    """_summary_

    Args:
        study_count (_type_): _description_
        study_n (_type_): _description_
        pop_count (_type_): _description_
        pop_n (_type_): _description_

    Returns:
        _type_: _description_
    """
    avar = study_count
    bvar = study_n - study_count
    cvar = pop_count - study_count
    dvar = pop_n - pop_count - bvar

    _, pval = stats.fisher_exact([[avar, bvar], [cvar, dvar]], alternative="greater")

    return pval


def binomialscipystats(study_count, study_n, pop_count, pop_n) -> float:
    k = study_count
    n = study_n
    p = pop_count / pop_n

    pval = stats.binomtest(k, n, p).pvalue

    return pval


def pvalue_calculate(study_count, study_n, pop_count, pop_n, method) -> float:
    if method == "fisher_scipy_stats":
        pval = fisherscipystats(study_count, study_n, pop_count, pop_n)
    elif method == "binomial_scipy_stats":
        pval = binomialscipystats(study_count, study_n, pop_count, pop_n)
    else:
        raise ValueError(f"{method} is not an available method for calculating pvalue")

    return pval

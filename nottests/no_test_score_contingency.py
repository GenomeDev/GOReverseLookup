"""
Scores basic contingency values
"""

from goreverselookup import fisher_exact_test

fisher_exact_test.compute_contingency(
    study_set=3,
    study_count=3,
    population_set=14,
    population_count=5
)

print("")

fisher_exact_test.compute_contingency(
    study_set=4,
    study_count=1,
    population_set=14,
    population_count=5
)
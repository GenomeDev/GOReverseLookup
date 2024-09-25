"""
Scores basic contingency values.

For GOReverseLookup: 
If you are testing whether a gene is more likely to be associated with a state of interest (SOI) than not, and you expect a higher occurrence in the study set, choose 'greater'.
If you are unsure of the direction but want to check for any difference in association, choose 'two-sided'.
If you expect the gene to occur less frequently in the study set, choose 'less'.

In GOReverseLookup, we are testing for genes that are overrepresented in the target SOIs,
thus we are choosing the 'greater' alternative hypothesis for Fisher's exact test.
"""

from goreverselookup import fisher_exact_test

fisher_exact_test.compute_contingency(
    study_set=3,
    study_count=3,
    population_set=14,
    population_count=5,
    hypothesis="greater"
)

print("")

fisher_exact_test.compute_contingency(
    study_set=4,
    study_count=1,
    population_set=14,
    population_count=5,
    hypothesis="greater"
)
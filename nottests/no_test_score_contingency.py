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

interest_set = 3
interest_count = 0
population_set = 28
population_count = 7

study_set = interest_set
study_count = interest_count
population_set = population_set
population_count = population_count

fisher_exact_test.compute_contingency(
    study_set=study_set,
    study_count=study_count,
    population_set=population_set,
    population_count=population_count,
    hypothesis="greater"
)

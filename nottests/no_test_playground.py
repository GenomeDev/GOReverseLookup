from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import GOTerm
from goreverselookup import GOApi
from goreverselookup import Cacher
from goreverselookup import JsonUtil

import time

# setup logger
import logging


# logger = logging.getLogger(__name__)
# logger.addHandler(logging.StreamHandler) # this doesn't work on windows
LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

"""
model = ReverseLookup.load_model("results/data.json")
nterms_score = nterms(model)
adv_prod_score = adv_product_score(model)
binom_score = binomial_test(model)
fisher_score = fisher_exact_test(model)
model.score_products(
    score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score]
)
model.model_settings.pvalue = 0.05  # set pvalue to be used in statistical analysis
model.perform_statistical_analysis(
    test_name="fisher_test", filepath="results/statistically_relevant_genes.json"
)
model.save_model("results/data.json")
"""


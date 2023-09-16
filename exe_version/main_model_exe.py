# This file demonstrates the same functionality as main_workflows_test.py, just using
# a ReverseLookup model. In the background, Workflows do the same function calls as you will
# see in this .py file.

import os
from goreverselookup import Cacher
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader

# setup logger
import logging

LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

logger.info("Starting Model Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

Cacher.init(cache_dir="cache")

# load the model from input file and query relevant data from the web
model = ReverseLookup.from_input_file("input_files/input.txt")
model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=0.1)
model.fetch_all_go_term_products(web_download=True, run_async=True, max_connections=35)
model.create_products_from_goterms()
model.fetch_ortholog_products(
    run_async=True, max_connections=10, semaphore_connections=4, req_delay=0.1
)
model.prune_products()
model.fetch_product_infos(
    refetch=False,
    run_async=True,
    max_connections=12,
    semaphore_connections=8,
    req_delay=0.1,
)
model.save_model("results/data.json")

# test model load from existing json, perform model scoring
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

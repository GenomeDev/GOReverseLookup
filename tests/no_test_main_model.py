# This file demonstrates the same functionality as main_workflows_test.py, just using
# a ReverseLookup model. In the background, Workflows do the same function calls as you will
# see in this .py file.

import os
import sys
from goreverselookup import Cacher, ModelStats
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import WebsiteParser

# setup logger
import logging

# logger = logging.getLogger(__name__)
# logger.addHandler(logging.StreamHandler)
LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

logger.info("Starting Model Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

Cacher.init(cache_dir="cache")
ModelStats.init()
WebsiteParser.init()

input_file = "input_files/input_rhartritis_noortho.txt"

# load the model from input file and query relevant data from the web
model = ReverseLookup.from_input_file(input_file)
# model.goterms = model.goterms[0:20]     # TODO: comment out
model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=1, max_connections=20)  # TODO: reenable this
model.fetch_all_go_term_products(web_download=True, run_async=False, delay=0.5, max_connections=5)
model.create_products_from_goterms()
model.products_perform_idmapping()
model.fetch_orthologs_products_batch_gOrth(target_taxon_number="9606")
model.fetch_ortholog_products(run_async=True, max_connections=20, semaphore_connections=10, req_delay=0.1)
model.prune_products()

#
# when using gorth_ortholog_fetch_for_indefinitive_orthologs as True,
# the ortholog count can go as high as 15.000 or even 20.000 -> fetch product infos
# disconnects from server, because we are seen as a bot.
# TODO: implement fetch_product_infos only for statistically relevant terms

# model.fetch_product_infos(
#    refetch=False,
#    run_async=True,
#    max_connections=15,
#    semaphore_connections=10,
#    req_delay=0.1,
# )

model.save_model("results/data.json")

# test model load from existing json, perform model scoring
model = ReverseLookup.load_model("results/data.json")
nterms_score = nterms(model)
adv_prod_score = adv_product_score(model)
binom_score = binomial_test(model)
fisher_score = fisher_exact_test(model)
model.score_products(score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score])

# model.model_settings.pvalue = 0.10  # set pvalue to be used in statistical analysis
model.perform_statistical_analysis(
	test_name="fisher_test", 
	filepath="results/statistically_relevant_genes.json", 
	exclude_opposite_regulation_direction_check=model.model_settings.exclude_opposite_regulation_direction_check
)
model.save_model("results/data.json")

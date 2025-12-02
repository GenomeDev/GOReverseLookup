# Rescores a model under a new p-value threshold

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

P_THRESHOLD = 0.00000005
MODEL_FILEPATH = "research_models\\BT_fish-growth\\direct_annotations\\0.05\\results\\data.json"
model = ReverseLookup.load_model(MODEL_FILEPATH)
model.model_settings.pvalue = P_THRESHOLD

GOTERM_NAME_FETCH_REQ_DELAY = 1.0
GOTERM_NAME_FETCH_MAX_CONNECTIONS = 20
GOTERM_GENE_FETCH_REQ_DELAY = 0.5   
GOTERM_GENE_FETCH_MAX_CONNECTIONS = 7
if 'goterm_name_fetch_req_delay' in model.model_settings.__dict__:
    GOTERM_NAME_FETCH_REQ_DELAY = model.model_settings.goterm_name_fetch_req_delay
if 'goterm_name_fetch_max_connections' in model.model_settings.__dict__:
    GOTERM_NAME_FETCH_MAX_CONNECTIONS = model.model_settings.goterm_name_fetch_max_connections
if 'goterm_gene_fetch_req_delay' in model.model_settings.__dict__:
    GOTERM_GENE_FETCH_REQ_DELAY = model.model_settings.goterm_gene_fetch_req_delay
if 'goterm_gene_fetch_max_connections' in model.model_settings.__dict__:
    GOTERM_GENE_FETCH_MAX_CONNECTIONS = model.model_settings.goterm_gene_fetch_max_connections

fisher_score = fisher_exact_test(model)
model.score_products(score_classes=[fisher_score])

# model.model_settings.pvalue = 0.10  # set pvalue to be used in statistical analysis
model.perform_statistical_analysis(
	test_name="fisher_test", 
	filepath="results/statistically_relevant_genes.json", 
)
model.save_model("results/data.json")

# This is the main file for the project that is used when goreverselookup is used
# from the command-line interface.

import os
import sys
from goreverselookup import Cacher, ModelStats
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import WebsiteParser

# setup logger
import logging
LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

logger.info("Starting GOReverseLookup analysis!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

def main(input_file:str):
    # Runs the GOReverseLookup analysis

    # setup
    Cacher.init(cache_dir="cache")
    ModelStats.init()
    WebsiteParser.init()
    
    # load the model from input file and query relevant data from the web
    model = ReverseLookup.from_input_file(input_file)
    model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=1, max_connections=20)  # TODO: reenable this
    model.fetch_all_go_term_products(web_download=True, run_async=True, delay=0.5, max_connections=20)
    model.fetch_all_go_term_products(web_download=True, run_async=True, delay=0.5, max_connections=20)
    model.create_products_from_goterms()
    model.products_perform_idmapping()
    model.fetch_orthologs_products_batch_gOrth(target_taxon_number="9606")
    model.fetch_ortholog_products(run_async=True, max_connections=20, semaphore_connections=10, req_delay=0.1)
    model.prune_products()
    model.save_model("results/data.json")

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

    # test model load from existing json, perform model scoring
    model = ReverseLookup.load_model("results/data.json")
    nterms_score = nterms(model)
    adv_prod_score = adv_product_score(model)
    binom_score = binomial_test(model)
    fisher_score = fisher_exact_test(model)
    model.score_products(score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score])
    model.perform_statistical_analysis(test_name="fisher_test", filepath="results/statistically_relevant_genes.json")
    # TODO: fetch info for stat relevant genes here
    model.save_model("results/data.json")

if len(sys.argv) != 2:
    print("Usage: goreverselookup <input_file>")
    sys.exit(1)
input_file = sys.argv[1]
logger.info(f"input_file = {input_file}")
main(input_file=input_file)



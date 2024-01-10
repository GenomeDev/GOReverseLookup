# This is the main file for the project that is used when goreverselookup is used from the command-line interface.

import argparse
import os
from goreverselookup import Cacher, ModelStats
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import WebsiteParser

# change directory to project root dir to include the logging_config.json file -> needed to setup logger
prev_cwd = os.getcwd()
project_root_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(project_root_dir)

# setup logger
import logging
LogConfigLoader.setup_logging_config(log_config_json_filepath="config/logging_config.json")
logger = logging.getLogger(__name__)

# change directory back to the main environment directory
os.chdir(prev_cwd)

logger.info("Starting GOReverseLookup analysis!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

def main(input_file:str, destination_dir:str = None):
    # Runs the GOReverseLookup analysis
    if destination_dir is None:
        destination_dir = os.path.dirname(os.path.realpath(input_file))

    # setup
    Cacher.init(cache_dir="cache")
    ModelStats.init()
    WebsiteParser.init()
    
    # load the model from input file and query relevant data from the web
    model = ReverseLookup.from_input_file(filepath=input_file, destination_dir=destination_dir)
    model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=1, max_connections=20)  # TODO: reenable this
    model.fetch_all_go_term_products(web_download=True, run_async=True, delay=0.5, max_connections=10)
    Cacher.save_data()
    model.create_products_from_goterms()
    model.products_perform_idmapping()
    Cacher.save_data()
    model.fetch_orthologs_products_batch_gOrth(target_taxon_number=f"{model.model_settings.target_organism.ncbi_id}") # TODO: change!
    model.fetch_ortholog_products(run_async=True, max_connections=15, semaphore_connections=7, req_delay=0.1)
    model.prune_products()
    model.bulk_ens_to_genename_mapping()
    model.save_model("results/data.json", use_dest_dir=True)

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
    model = ReverseLookup.load_model("results/data.json", destination_dir=destination_dir)
    nterms_score = nterms(model)
    adv_prod_score = adv_product_score(model)
    binom_score = binomial_test(model)
    fisher_score = fisher_exact_test(model)
    model.score_products(score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score])
    model.perform_statistical_analysis(test_name="fisher_test", filepath="results/statistically_relevant_genes.json", use_dest_dir=True)
    # TODO: fetch info for stat relevant genes here
    model.save_model("results/data.json", use_dest_dir=True)

#if len(sys.argv) != 2:
#    print("Usage: goreverselookup <input_file>")
#    sys.exit(1)
#input_file = sys.argv[1]
#logger.info(f"input_file = {input_file}")

parser = argparse.ArgumentParser(description="Usage: goreverselookup <input_file_path> --<destination_directory> ('--' denotes an optional parameter)")
parser.add_argument('input_file', help="The absolute path to the input file for GOReverseLookup analysis.")
parser.add_argument('--destination_dir', help="The directory where output and intermediate files will be saved. If unspecified, output directory will be selected as the root directory of the supplied input file.")
# TODO: debug arguments

# parse the command-line arguments
args = parser.parse_args()
input_file = args.input_file
destination_dir = args.destination_dir
main(input_file=input_file, destination_dir=destination_dir)



from goreverselookup import ReverseLookup
from goreverselookup import WorkflowTwo
from goreverselookup import Cacher
from goreverselookup import JsonUtil

from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import GOAnnotationsFile

import logging
from logging import config
log_config_json_filepath = "app/goreverselookup/src/logging_config.json"
log_config_dict = JsonUtil.load_json(log_config_json_filepath)
config.dictConfig(log_config_dict)
logger = logging.getLogger(__name__)

Cacher.init(cache_dir="app/goreverselookup/cache")

# *** Workflows algorithm ***
workflow = WorkflowTwo(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()

# *** Model algorithm ***
#model = ReverseLookup.from_input_file("input_files/input.txt")
#model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=0.1)
#model.fetch_all_go_term_products(web_download=True, run_async=True)
#model.create_products_from_goterms()
#model.fetch_ortholog_products(run_async=True, max_connections=15, semaphore_connections=5, req_delay=0.1)
#model.prune_products()
#model.save_model("results/data.json")

# perform model scoring and also test model load from json functionality
#model = ReverseLookup.load_model("results/data.json")
#nterms_score = nterms(model)
#adv_prod_score = adv_product_score(model)
#binom_score = binomial_test(model)
#fisher_score = fisher_exact_test(model)
#model.score_products(score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score])
#model.model_settings.pvalue = 0.05 # set pvalue to be used in statistical analysis
#model.perform_statistical_analysis(test_name="fisher_test", filepath="results/statistically_relevant_genes.json")
#model.save_model("results/data.json")

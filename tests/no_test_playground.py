from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import GOTerm
from goreverselookup import GOApi

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

model_old = ReverseLookup.load_model("results/run2/data.json")
model_new = ReverseLookup.load_model("results/data.json")

# model_old = source, model_new = reference
mismatches = model_old.compare_to(model_new, compare_field="goterms", compare_subfields=["products"])

file = open(f"log_output/goterm_product_mismatches.txt", "a+")

for mismatch in mismatches:
    # mismatch is a GO id
    if "No source element" in mismatches[mismatch]["mismatches"][0]:
        logger.info(f"Source GO Term {mismatch} was not found in model_new!")
        file.write(f"Source GO Term {mismatch} was not found in model_new!\n")
        continue
    
    split = mismatches[mismatch]['mismatches'][0].split("\\n")
      
    missing_reference_elements_in_source = split[1].split("in src: ")[1] 
    missing_source_elements_in_reference = split[2].split("in reference: ")[1]
    # convert from string back to list
    missing_reference_elements_in_source = eval(missing_reference_elements_in_source)
    missing_source_elements_in_reference = eval(missing_source_elements_in_reference)

    products_model_old = model_old.get_goterm(mismatch).products
    products_model_new = model_new.get_goterm(mismatch).products

    logger.info(f"{mismatch}: {len(products_model_old)} products in model_old, {len(products_model_new)} products in model_new")
    logger.info(f"  - products present in model_old, but not in model_new ({len(missing_source_elements_in_reference)}): {missing_source_elements_in_reference}")
    logger.info(f"  - products present in model_new, but not in model_old ({len(missing_reference_elements_in_source)}): {missing_reference_elements_in_source}")
    

    file.write(f"{mismatch}: {len(products_model_old)} products in model_old, {len(products_model_new)} products in model_new\n")
    file.write(f"  - products present in model_old, but not in model_new ({len(missing_source_elements_in_reference)}): {missing_source_elements_in_reference}\n")
    file.write(f"  - products present in model_new, but not in model_old ({len(missing_reference_elements_in_source)}): {missing_reference_elements_in_source}\n")
    file.write(f"\n")
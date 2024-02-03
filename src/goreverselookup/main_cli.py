# This is the main file for the project that is used when goreverselookup is used from the command-line interface.

import argparse
import os
import sys
from goreverselookup import Cacher, ModelStats
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import WebsiteParser
from goreverselookup import JsonUtil, FileUtil
import pandas as pd

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

def generate_report(results_file, model_data):
    # TODO: more refined reporting functionality

    if isinstance(results_file, str):
        results_file = JsonUtil.load_json(results_file)
    if isinstance(model_data, str):
        model_data = JsonUtil.load_json(model_data)
        
    print(f"p-value: {model_data['model_settings']['pvalue']}")
    print(f"include indirect annotations: {model_data['model_settings']['include_indirect_annotations']}")
    # TODO: evidence codes
    # TODO: other settings

    target_SOIs = model_data['target_SOIs']
    SOIs = [] # e.g. ['chronic_inflammation+', 'cancer+'] # todo: also account for reverse tSOIs (eg negative)
    stat_rev_key = ""
    if len(target_SOIs) == 1:
        SOIs.append(f"{target_SOIs[0]['SOI']}+")
        SOIs.append(f"{target_SOIs[0]['SOI']}-")
        stat_rev_key=f"{target_SOIs[0]['SOI']}{target_SOIs[0]['direction']}"
    else:
        for tSOI in target_SOIs:
            intermediate_stat_rev_key = f"{tSOI['SOI']}{tSOI['direction']}"
            SOIs.append(f"{tSOI['SOI']}+")
            SOIs.append(f"{tSOI['SOI']}-")
            stat_rev_key = intermediate_stat_rev_key if stat_rev_key == "" else f"{stat_rev_key}:{intermediate_stat_rev_key}"
        
    stat_rev_genes = results_file[stat_rev_key]
        
    #separator = " "
    separator = "\t"
        
    pdata = { # data for pandas export to excel
        'genename': []
    }

    # table start text
    stext = f"gene{separator}"
    for SOI in SOIs:
        stext = f"{stext}{separator}{SOI}" 
        pdata[SOI] = [] # init SOI data for pandas export to excel
    print(stext)
        
    # generate table between genes and tr_SOIs
    for sr_gene in stat_rev_genes:
        genename = sr_gene['genename']
        pdata['genename'].append(genename)
        pvalues = []
        for SOI in SOIs:
            if SOI in sr_gene['scores']['fisher_test']:
                if 'pvalue_corr' in sr_gene['scores']['fisher_test'][SOI]:
                    pvalue = sr_gene['scores']['fisher_test'][SOI]['pvalue_corr']
                    pvalue_form = "{:.2e}".format(pvalue)
                else:
                    pvalue_form = "/"
            else:
                pvalue_form = "/"
            pvalues.append(pvalue_form)
            pdata[SOI].append(pvalue_form)
        rowtext = f"{genename}{separator}"
        for pval in pvalues:
            rowtext = f"{rowtext}{separator}{pval}"
        print(rowtext)
        
    # to excel
    root = FileUtil.backtrace(input_file, 1)
    stat_rev_genes_xlsx_path = f"{root}/stat_rev_genes.xlsx"
    df = pd.DataFrame(pdata)
    df.to_excel(stat_rev_genes_xlsx_path, index=True, header=True)

    sys.exit(0)

def main(input_file:str, destination_dir:str = None, report:bool = False, model_data_filepath:str = None):
    logger.info(f"Starting GOReverseLookup analysis with input params:")
    logger.info(f"  - input_file: {input_file}")
    logger.info(f"  - destination_dir: {destination_dir}")
    logger.info(f"  - report: {report}")
    logger.info(f"  - model_data_filepath: {model_data_filepath}")

    model_data = None
    if model_data_filepath is not None:
        model_data = JsonUtil.load_json(model_data_filepath)
    
    if report is True and model_data is not None: # should generate report only
        results_file = JsonUtil.load_json(input_file)
        generate_report(results_file, model_data)
         
    # Runs the GOReverseLookup analysis
    if destination_dir is None:
        destination_dir = os.path.dirname(os.path.realpath(input_file))

    # setup
    Cacher.init(cache_dir="cache")
    ModelStats.init()
    WebsiteParser.init()
    
    # load the model from input file and query relevant data from the web
    if model_data is None:
        model = ReverseLookup.from_input_file(filepath=input_file, destination_dir=destination_dir)
        print(f"Model was created from input file: {input_file}")
    else:
        # model_data[''] # TODO ADD DESTINATION DIR HERE !!!!!
        model = ReverseLookup.load_model(model_data_filepath, destination_dir=destination_dir)
        print(f"Model was created from a previous model_data dictionary: {model_data_filepath}")
    model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=1, max_connections=20) 
    model.fetch_all_go_term_products(web_download=True, run_async=True, delay=0.5, max_connections=7)
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

    # TODO
    # generate_report("results/statistically_relevant_genes.json", "results/data.json")

#if len(sys.argv) != 2:
#    print("Usage: goreverselookup <input_file>")
#    sys.exit(1)
#input_file = sys.argv[1]
#logger.info(f"input_file = {input_file}")

# UNCOMMENT THIS SECTION FOR PRODUCTION CODE !!!

parser = argparse.ArgumentParser(description="Usage: goreverselookup <input_file_path> --<destination_directory> ('--' denotes an optional parameter)")
parser.add_argument('input_file', help="The absolute path to the input file for GOReverseLookup analysis or to the resulting file if used with the --report optional parameter.")
parser.add_argument('--destination_dir', help="The directory where output and intermediate files will be saved. If unspecified, output directory will be selected as the root directory of the supplied input file.")
parser.add_argument('--report', help="Values: True or False. Specify this optional parameter to generate a report of statistically significant genes (the input file must point to a statistically_significant_genes.json)")
parser.add_argument('--model_datafile', help="The main research model data file path (usually generated as data.json). If specifying model_datafile, it will create the research model from the supplied model datafile (precedence over the input file). If left unspecified and using '--report True', then an attempt is made to infer model_datafile from the root directory of input_filepath. Thus, if statistically_significant_genes.json and data.json are saved in the same directory, --report True can be ran without the model_datafile parameter.")
# TODO: debug arguments

# parse the command-line arguments
args = parser.parse_args()
input_file = args.input_file
destination_dir = args.destination_dir

report = args.report
if report is not None and report.upper() == "TRUE":
    report = True
else:
    report = False
    
model_data_filepath = args.model_datafile
if model_data_filepath is None:
    print("No model data filepath was specified, auto-inferring model data.")
    root = FileUtil.backtrace(input_file, 1) # move 1 file up to root dir
    m_data_filepath = os.path.join(root, "data.json")
    m_data_filepath = m_data_filepath.replace("\\", "/")
    if FileUtil.check_path(m_data_filepath, auto_create=False):
        if FileUtil.get_file_size(m_data_filepath, "kb") > 5: # if ReverseLookup data file is greater than 5kb then assign, otherwise it's most likely an error
            print(f"Model data filepath found: {m_data_filepath}")
            model_data_filepath = m_data_filepath
    else:
        print(f"Model data not found. Attempted file search at {model_data_filepath}")


# test arguments for debugging, remove these
#input_file = "F:\\Development\\python_environments\\goreverselookup\\research_models\\chr_infl_cancer\\ind_ann,p=5e-8,IEA+\\input.txt"
#destination_dir = None
#report = False
#model_data_filepath = None

main(input_file=input_file, destination_dir=destination_dir, report=report, model_data_filepath=model_data_filepath)



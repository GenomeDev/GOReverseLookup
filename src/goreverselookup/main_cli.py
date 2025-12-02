import argparse
import os
import copy
import logging
from logging.handlers import RotatingFileHandler

from goreverselookup import Cacher, ModelStats
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import WebsiteParser
from goreverselookup import JsonUtil, FileUtil
import pandas as pd


# ---------------------------
# Logging configuration
# ---------------------------

def _reset_root_logger():
    """Remove all handlers from the root logger to avoid duplicate logs when reconfiguring."""
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass

def configure_logging_for_input_file(input_file: str, to_console: bool = True, level: int = logging.DEBUG):
    """
    Configure the *root* logger so that ALL loggers from every module (logging.getLogger(__name__))
    propagate here and are written into log.log located in the *root folder of the input file*.

    This replaces any prior configuration to ensure clean per-input logging.
    """
    # Resolve root folder for the input file
    root_dir = os.path.dirname(os.path.realpath(input_file))
    log_path = os.path.join(root_dir, "log.log")

    # Ensure directory exists (it should) and that we can write the file
    os.makedirs(root_dir, exist_ok=True)

    # Reset any existing handlers to prevent duplicate entries
    _reset_root_logger()

    # Build handlers
    file_handler = RotatingFileHandler(
        filename=log_path,
        mode="a",
        maxBytes=10 * 1024 * 1024,  # 10 MB
        backupCount=3,
        encoding="utf-8",
        delay=False,
    )

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler.setFormatter(formatter)
    file_handler.setLevel(level)

    handlers = [file_handler]

    if to_console:
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        console.setLevel(logging.INFO)
        handlers.append(console)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    for h in handlers:
        root_logger.addHandler(h)

    # Make sure warnings captured too
    logging.captureWarnings(True)

    # Slightly quiet down very chatty third-party libs if needed:
    logging.getLogger("asyncio").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)

    # Return a module-specific logger for this file
    return logging.getLogger(__name__)


logger = logging.getLogger(__name__)  # Will be re-bound after parsing args by configure_logging_for_input_file


# ---------------------------
# Core functionality
# ---------------------------

def generate_report(results_file: str, model_data):
    # more refined reporting functionality
    logger.info("Generating report...")
    results = JsonUtil.load_json(results_file)
    if isinstance(model_data, str):
        model_data = JsonUtil.load_json(model_data)

    logger.info("Model settings:")
    logger.info(f"  - p-value: {model_data['model_settings']['pvalue']}")
    logger.info(f"  - include indirect annotations: {model_data['model_settings']['include_indirect_annotations']}")

    target_SOIs = model_data['target_SOIs']
    SOIs = []  # e.g. ['chronic_inflammation+', 'cancer+']
    stat_rev_key = ""
    if len(target_SOIs) == 1:
        SOIs.append(f"{target_SOIs[0]['SOI']}+")
        SOIs.append(f"{target_SOIs[0]['SOI']}-")
        stat_rev_key = f"{target_SOIs[0]['SOI']}{target_SOIs[0]['direction']}"
    else:
        for tSOI in target_SOIs:
            intermediate_stat_rev_key = f"{tSOI['SOI']}{tSOI['direction']}"
            SOIs.append(f"{tSOI['SOI']}+")
            SOIs.append(f"{tSOI['SOI']}-")
            stat_rev_key = intermediate_stat_rev_key if stat_rev_key == "" else f"{stat_rev_key}:{intermediate_stat_rev_key}"

    stat_rev_genes = results[stat_rev_key]

    separator = "\t"

    pdata = {  # data for pandas export to excel
        'genename': [],
        'id_synonym': []
    }

    # table header text
    header = f"gene{separator}" + separator.join(SOIs)
    logger.info(header)

    # generate table between genes and tr_SOIs
    for sr_gene in stat_rev_genes:
        genename = sr_gene['genename']
        pdata['genename'].append(genename)
    
        id_synonyms = sr_gene.get('id_synonyms', [])
        if isinstance(id_synonyms, list) and len(id_synonyms) > 0:
            first_id_synonym = id_synonyms[0]
        elif isinstance(id_synonyms, str):
            # in case it is already a single string for some reason
            first_id_synonym = id_synonyms
        else:
            first_id_synonym = "/"
        pdata['id_synonym'].append(first_id_synonym)
        
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
            pdata[SOI] = pdata.get(SOI, [])
            pdata[SOI].append(pvalue_form)
        rowtext = f"{genename}{separator}" + separator.join(pvalues)
        logger.info(rowtext)

    # to excel
    root = FileUtil.backtrace(results_file, 1)
    stat_rev_genes_xlsx_path = f"{root}/stat_rev_genes.xlsx"
    df = pd.DataFrame(pdata)
    df.to_excel(stat_rev_genes_xlsx_path, index=True, header=True)
    logger.info(f"Results saved to: {stat_rev_genes_xlsx_path}")


def main(input_file: str, destination_dir: str = None, report: bool = False, model_data_filepath: str = None,
         rescore_model: ReverseLookup = None):
    logger.info("Starting GOReverseLookup analysis with input params:")
    logger.info(f"  - input_file: {input_file}")
    logger.info(f"  - destination_dir: {destination_dir}")
    logger.info(f"  - report: {report}")
    logger.info(f"  - model_data_filepath: {model_data_filepath}")
    logger.debug(f"  - os.getcwd() = {os.getcwd()}")

    model_data = None
    if model_data_filepath is not None:
        model_data = JsonUtil.load_json(model_data_filepath)
    else:
        logger.info(f"Model data filepath was None. Attempting auto infer data.json from {input_file}")
        root = FileUtil.backtrace(input_file, 1)  # move 1 file up to root dir 
        m_data_filepath = os.path.join(root, "data.json").replace("\\", "/")
        logger.debug(f"m_data_filepath={m_data_filepath}")
        if FileUtil.check_path(m_data_filepath, auto_create=False):
            if FileUtil.get_file_size(m_data_filepath, "kb") > 0:
                logger.info(f"Model data filepath was found by auto infer: {m_data_filepath}")
                model_data_filepath = m_data_filepath
                model_data = JsonUtil.load_json(model_data_filepath)
        else:
            logger.warning("Model data was not found by auto-infer.")

    if report is True and model_data is not None:  # should generate report only
        generate_report(results_file=input_file, model_data=model_data)
        return
    elif report is True and model_data is None:  # error
        logger.error(
            "Report is True, but no model data (data.json) files were found. "
            "You need to keep both statistically_relevant_genes.json and data.json in the same folders without deleting them!"
        )
        return

    if destination_dir is None:
        destination_dir = os.path.dirname(os.path.realpath(input_file))

    # if True, only perform rescoring with new values
    if rescore_model is not None:
        logger.info("Rescoring mode enabled — updating p-value/indirect annotations only.")
        new_model = ReverseLookup.from_input_file(filepath=input_file,
                                                  destination_dir=destination_dir) if model_data is None else ReverseLookup.load_model(
            model_data_filepath, destination_dir=destination_dir)
        rescore_model_copy = copy.deepcopy(rescore_model)
        rescore_model_copy.model_settings = new_model.model_settings
        rescore_model_copy.statistically_relevant_products = {}  # reset
        fisher_score = fisher_exact_test(rescore_model_copy)
        rescore_model_copy.score_products(score_classes=[fisher_score])
        rescore_model_copy.perform_statistical_analysis(
            test_name="fisher_test",
            filepath="results/statistically_relevant_genes.json",
            use_dest_dir=True,
            two_tailed=rescore_model.model_settings.two_tailed
        )
        rescore_model_copy.save_model("results/data.json", use_dest_dir=True)
        return

    # setup
    Cacher.init(cache_dir="cache")
    ModelStats.init()
    WebsiteParser.init()

    # load the model from input file and query relevant data from the web
    if model_data is None:
        model = ReverseLookup.from_input_file(filepath=input_file, destination_dir=destination_dir)
        logger.info(f"Model was created from input file: {input_file}")
    else:
        model = ReverseLookup.load_model(model_data_filepath, destination_dir=destination_dir)
        logger.info(f"Model was created from a previous model_data dictionary: {model_data_filepath}")

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

    model.fetch_all_go_term_names_descriptions(run_async=model.model_settings.goterm_name_fetch_async, req_delay=GOTERM_NAME_FETCH_REQ_DELAY,
                                               max_connections=GOTERM_NAME_FETCH_MAX_CONNECTIONS)
    model.fetch_all_go_term_products(web_download=True, run_async=model.model_settings.goterm_gene_fetch_async, delay=GOTERM_GENE_FETCH_REQ_DELAY,
                                     max_connections=GOTERM_GENE_FETCH_MAX_CONNECTIONS,
                                     request_params={"rows": 10000000})
    Cacher.save_data()
    model.create_products_from_goterms()
    model.products_perform_idmapping()
    model.fetch_orthologs_products_batch_gOrth(target_taxon_number=f"{model.model_settings.target_organism.ncbi_id}")  # TODO: change!
    model.fetch_ortholog_products(run_async=True, max_connections=15, semaphore_connections=7, req_delay=0.1)
    model.prune_products()
    model.bulk_ens_to_genename_mapping()
    model.save_model("results/data.json", use_dest_dir=True)
    Cacher.save_data()

    # Perform scoring
    model = ReverseLookup.load_model("results/data.json", destination_dir=destination_dir)
    fisher_score = fisher_exact_test(model)
    model.score_products(score_classes=[fisher_score])
    model.perform_statistical_analysis(
        test_name="fisher_test",
        filepath="results/statistically_relevant_genes.json",
        use_dest_dir=True,
        two_tailed=model.model_settings.two_tailed
    )
    model.save_model("results/data.json", use_dest_dir=True)
    return model

# ---------------------------
# CLI
# ---------------------------

parser = argparse.ArgumentParser(
    description="Usage: goreverselookup <input_file_path> --<destination_directory> ('--' denotes an optional parameter)"
)
parser.add_argument('input_file',
                    help="The absolute path to the input file for GOReverseLookup analysis or to the resulting file if used with the --report optional parameter.")
parser.add_argument('--destination_dir',
                    help="The directory where output and intermediate files will be saved. If unspecified, output directory will be selected as the root directory of the supplied input file.")
parser.add_argument('--report',
                    help="Values: True or False. Specify this optional parameter to generate a report of statistically significant genes (the input file must point to a statistically_significant_genes.json)")
parser.add_argument('--model_datafile',
                    help="The main research model data file path (usually generated as data.json). If specifying model_datafile, it will create the research model from the supplied model datafile (precedence over the input file). If left unspecified and using '--report True', then an attempt is made to infer model_datafile from the root directory of input_filepath. Thus, if statistically_significant_genes.json and data.json are saved in the same directory, --report True can be ran without the model_datafile parameter.")
parser.add_argument("--full_directory_op",
                    help="Specify the root directory, all subdirectories will be scanned for a file named as 'input_file' and the operation will be performed on all these files.")
parser.add_argument("--rescore",
                    help="Used in conjunction with full_directory_op. If True, will preserve the initially computed research model, but will only update the p-value, indirect annotations. WARNING: if using this setting with models with different evidence code settings, there will be errors!")
# TODO: debug arguments

# parse the command-line arguments
args = parser.parse_args()
input_file = args.input_file
destination_dir = args.destination_dir

report = False
rescore = False
if args.report is not None:
    report = True if args.report.upper() == "TRUE" else False
if args.rescore is not None:
    rescore = True if args.rescore.upper() == "TRUE" else False
static_rescore_model = None  # model used for re-scoring

full_directory_op = args.full_directory_op

model_data_filepath = args.model_datafile
if model_data_filepath is None:
    # NOTE: Keep this discovery via logging (not print) so it also lands in log.log after configuration.
    pass

# Build list of input files (single or directory sweep)
input_files = []
if full_directory_op is None:
    input_files = [input_file]
else:
    # infer all input files from the specified input file and the root directory provided in full_directory_op
    input_file_name = FileUtil.get_filename(input_file)
    input_files = FileUtil.get_directory_files(full_directory_op, input_file_name)

# Iterate files; reconfigure logging for each file so each root gets its own log.log
for idx, f in enumerate(input_files):
    # Configure logging per input file root
    logger = configure_logging_for_input_file(f)

    logger.info("==============================================")
    logger.info(f"Processing file [{idx+1}/{len(input_files)}]: {f}")
    logger.info("==============================================")

    # If model_data_filepath not provided, attempt auto-infer (now logged)
    if model_data_filepath is None:
        logger.info("No model data filepath was specified, auto-inferring model data.")
        root = FileUtil.backtrace(f, 1)  # move 1 file up to root dir
        m_data_filepath = os.path.join(root, "data.json").replace("\\", "/")
        if FileUtil.check_path(m_data_filepath, auto_create=False):
            if FileUtil.get_file_size(m_data_filepath, "kb") > 5:
                logger.info(f"Model data filepath found: {m_data_filepath}")
                model_data_fp_for_this_run = m_data_filepath
            else:
                logger.warning(f"Found data.json but size too small to trust at: {m_data_filepath}")
                model_data_fp_for_this_run = None
        else:
            logger.warning(f"Model data not found. Attempted file search at {m_data_filepath}")
            model_data_fp_for_this_run = None
    else:
        model_data_fp_for_this_run = model_data_filepath

    model = main(
        input_file=f,
        destination_dir=destination_dir,
        report=report,
        model_data_filepath=model_data_fp_for_this_run,
        rescore_model=static_rescore_model
    )

    # if destination_dir is none, make it the directory of the input file (f)
    if destination_dir is None:
        destination_dir = os.path.dirname(f)

    dest_results = os.path.join(destination_dir, "results", "statistically_relevant_genes.json")
    dest_data = os.path.join(destination_dir, "results", "data.json")
    generate_report(results_file=dest_results, model_data=dest_data)

    if rescore and idx == 0:
        static_rescore_model = model

logger.info("All tasks completed.")
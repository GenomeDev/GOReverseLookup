# This file demonstrates how to run the analysis using Workflows.

# import
import os
from goreverselookup import PrimaryWorkflow
from goreverselookup import Cacher
from goreverselookup import LogConfigLoader
# setup logger
import logging

# logger = logging.getLogger(__name__)
# logger.addHandler(logging.StreamHandler) # this doesn't work on windows
LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

logger.info("Starting Workflows Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

# setup and run workflows
Cacher.init(cache_dir="cache")
# Cacher.clear_cache("ALL")
workflow = PrimaryWorkflow(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()

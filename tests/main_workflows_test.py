# This file demonstrates how to run the analysis using Workflows.

import os
from goreverselookup import PrimaryWorkflow
from goreverselookup import Cacher
from goreverselookup import LogConfigLoader

# setup logger
import logging
LogConfigLoader.setup_logging_config()
logger = logging.getLogger(__name__)

logger.info(f"Starting Workflows Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

# setup and run workflows
Cacher.init(cache_dir="app/goreverselookup/cache")
# Cacher.clear_cache("ALL")
workflow = PrimaryWorkflow(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()
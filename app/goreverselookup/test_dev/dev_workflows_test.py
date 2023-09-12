# This file demonstrates how to run the analysis using Workflows.

import os
from ..src.Workflows import PrimaryWorkflow
from ..src.util.CacheUtil import Cacher
from ..src.util.LogConfigLoader import LogConfigLoader

# setup logger
import logging
LogConfigLoader.setup_logging_config()
logger = logging.getLogger(__name__)

logger.info(f"Starting Workflows Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

# setup and run workflows
Cacher.init(cache_dir="app/goreverselookup/cache")
workflow = PrimaryWorkflow(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()
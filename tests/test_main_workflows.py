# This file demonstrates how to run the analysis using Workflows.

import os
from goreverselookup import PrimaryWorkflow
from goreverselookup import Cacher

# setup logger
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler)

logger.info(f"Starting Workflows Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

# setup and run workflows
Cacher.init(cache_dir="cache")
# Cacher.clear_cache("ALL")
workflow = PrimaryWorkflow(input_file_fpath="tests/test_input.txt", save_folder_dir="results")
workflow.run_workflow()
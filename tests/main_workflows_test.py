# This file demonstrates how to run the analysis using Workflows.

import os
from goreverselookup import WorkflowTwo
from goreverselookup import Cacher
from goreverselookup import JsonUtil

# setup logger
import logging
from logging import config
log_config_json_filepath = "app/goreverselookup/src/logging_config.json"
log_config_dict = JsonUtil.load_json(log_config_json_filepath)
config.dictConfig(log_config_dict)
logger = logging.getLogger(__name__)

logger.info(f"Starting Workflows Test!")
logger.info(f"os.getcwd() =  {os.getcwd()}")

# setup and run workflows
Cacher.init(cache_dir="app/goreverselookup/cache")
workflow = WorkflowTwo(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()
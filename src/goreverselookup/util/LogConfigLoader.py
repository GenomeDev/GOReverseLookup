# Previously, to setup the logger in a class with a config file, you needed to
# perform the following calls:
#
#   from goreverselookup import JsonUtil
#   import logging
#   from logging import config
#
#   log_config_json_filepath = "app/goreverselookup/src/logging_config.json"
#   log_config_dict = JsonUtil.load_json(log_config_json_filepath)
#   config.dictConfig(log_config_dict)
#   logger = logging.getLogger(__name__)
#
# Using LogConfigLoader, you reduce the complexity like so:
#
#   from goreverselookup import LogConfigLoader
#   import logging
#
#   LogConfigLoader.setup_logging_config()
#   logger = logging.getLogger(__name__)


from .JsonUtil import JsonUtil
from logging import config
import logging
logger = logging.getLogger(__name__)
# from goreverselookup import logger


class LogConfigLoader:
    def __init__():
        """ """
        pass

    @classmethod
    def setup_logging_config(
        cls, log_config_json_filepath="data_files/logging_config.json"
    ):
        """
        Performs the following function calls in order to setup the logging configuration:
            log_config_dict = JsonUtil.load_json(log_config_json_filepath)
            config.dictConfig(log_config_dict)
        """
        log_config_dict = JsonUtil.load_json(log_config_json_filepath, display_json=True)
        config.dictConfig(log_config_dict)
        logger.info(f"Setup log config using: {log_config_json_filepath}")


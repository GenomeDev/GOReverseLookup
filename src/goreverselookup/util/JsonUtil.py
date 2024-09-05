import json
import os
from types import SimpleNamespace
import platform

from .FileUtil import FileUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class JsonUtil:
    @classmethod
    def load_json(cls, filepath: str, display_json:bool=False):
        """
        Loads a json file and returns the json object (a dictionary).

        Parameters:
          - (str) filepath
        """
        logger.info(f"Load JSON received input filepath: {filepath}")
        initial_filepath = filepath

        if not os.path.exists(filepath):
            FileUtil.create_empty_file(filepath)
            
            if not os.path.isabs(filepath) and not platform.system == 'Darwin': # bugfix: fileutil.find_file doesn't work on Mac (platform.system == 'Darwin')
                fileutil = FileUtil()
                filepath = fileutil.find_file(filepath)  # attempt backtrace file search
                logger.info(f"Filepath after file search: {filepath}")

                if filepath is None:
                    raise Exception(f"Filepath when attempting load JSON is None! Initial filepath was {initial_filepath}")
                
        # bugfix: if filepath is empty, I want load_json to return {} instead of JSONDecodeError
        if FileUtil.is_file_empty(filepath):
            return {}

        with open(filepath, "r") as f:
            data = json.load(f)
            if display_json:
                logger.info(f"JSON visualisation:")
                logger.info(data)
            return data

    @classmethod
    def save_json(cls, data_dictionary: dict, filepath: str):
        """
        Saves the data_dictionary as json to the filepath.

        Parameters:
          - (dict) data_dictionary: the data to be saved as the json
          - (str) filepath: the filepath where the json is to be stored
        """
        if ".json" not in filepath:
            filepath = f"{filepath}.json"
        
        FileUtil.check_path(filepath)

        logger.info(f"Saving json to: {filepath}")
        with open(filepath, "w") as f:
            json.dump(data_dictionary, f, indent=4)

    @classmethod
    def class_to_json(cls, var):
        json_data = {}
        for attr_name, attr_value in vars(var).items():
            if attr_value is None:
                # this happens for example when only attr_name is defined without a value; e.g. "ortholog_organisms" without any following ortholog organisms.
                continue
            # custom handling for target_organism and ortholog_organisms, as they are code objects -> convert them to json
            if not callable(attr_value) and not attr_name.startswith("__"):
                # append to json_data result dict
                json_data[attr_name] = attr_value
        return json_data


class JsonToClass:
    object_representation = ""
    source_json = ""

    def __init__(self, data):
        """
        Converts a JSON eg. '{"name": "John Smith", "hometown": {"name": "New York", "id": 123}}'
        to a class, which can be queried with eg. x.name, x.hometown.name, x.hometown.id

        Warning: keys and values must be double-quote delimited (single-quotes are bugged)
        Warning: avoid usage of single quotes inside json string values (doing so breaks this code)

        Tip: http://json.parser.online.fr/ is very useful for debugging JSON errors, just remove starting ' and '

        Params:
            - data: a json representation of the data

        Returns: a class with the following member fields
            - object_representation: a SimpleNamespace instance, representing the objectified json
            - source_json: the json from which the SimpleNamespace was built

        Example usage: you want to convert a JSON to a Python object (class instance):
        json_string = '{"name": "John Smith", "hometown": {"name": "New York", "id": 123}}'
        json_namespace = SimpleNamespaceCustom(json_string)
        object_representation = json_namespace.object_representation
        src_json = json_namespace.source_json
        """

        # Parse JSON into an object with attributes corresponding to dict keys.
        if (
            "isn't" in data
        ):  # TODO: this is a hardcoded solution. Avoid any ' characters in the json values or change json to class loading.
            data = data.replace("isn't", "is not")

        if "None" in data:  # hardcoded bugfix
            data = data.replace("None", '"None"')

        if "inf" in data:  # hardcoded bugfix
            data = data.replace("inf", '"inf"')

        if "'" in data:
            data = data.replace("'", '"')  # SimpleNamespace bugs out with single-quotes

        if "nan" in data:
            data = data.replace("nan", "0")

        x = json.loads(data, object_hook=lambda d: SimpleNamespace(**d))
        self.object_representation = x
        self.source_json = data


class SimpleNamespaceUtil:
    """
    Utility functions for SimpleNamespace
    """

    def __init__():
        return 0

    @staticmethod
    def simpleNamespace_to_json(simple_namespace: SimpleNamespace):
        """
        TODO: incomplete function
        Converts a simpleNamespace object to a json string
        """
        # def iterate_members(base_member, result={}):
        #    submembers = []
        #    if isinstance(base_member, object):
        #        for attr in dir(base_member):
        #            if not callable(getattr(base_member, attr)) and not attr.startswith("__"):
        #                value = getattr(base_member, attr)
        #                if isinstance(value, object):
        #                    iterate_members(value, result)
        #                else:
        #                    result[attr] = value
        #                submembers.append(attr)
        #
        # member_fields = []
        # for attr in dir(simple_namespace):
        #    if not callable(getattr(simple_namespace, attr)) and not attr.startswith("__"):
        #        member_fields.append(attr)
        #
        return 0


def json_to_class(data: str):
    """
    Converts a JSON eg. '{"name": "John Smith", "hometown": {"name": "New York", "id": 123}}'
    to a class, which can be queried with eg. x.name, x.hometown.name, x.hometown.id

    Warning: keys and values must be double-quote delimited (single-quotes are bugged)
    Warning: avoid usage of single quotes inside json string values (doing so breaks this code)

    Tip: http://json.parser.online.fr/ is very useful for debugging JSON errors, just remove starting ' and '

    Params:
      - data: a json representation of the data

    Returns:
      - a class from the input json data
    """
    if (
        "isn't" in data
    ):  # TODO: this is a hardcoded solution. Avoid any ' characters in the json values or change json to class loading.
        data = data.replace("isn't", "is not")

    if "None" in data:  # hardcoded bugfix
        data = data.replace("None", '"None"')

    if "inf" in data:  # hardcoded bugfix
        data = data.replace("inf", '"inf"')

    if "'" in data:
        data = data.replace("'", '"')  # SimpleNamespace bugs out with single-quotes

    if "nan" in data:
        data = data.replace("nan", "0")

    x = json.loads(data, object_hook=lambda d: SimpleNamespace(**d))
    return x

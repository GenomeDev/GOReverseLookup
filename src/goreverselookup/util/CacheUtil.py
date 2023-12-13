import atexit
import os

from .JsonUtil import JsonUtil
from .FileUtil import FileUtil
from .Timer import Timer

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class Cacher:
    CACHE_FILEPATH_URLS = ""# filepath to the file containing online url queries and the URL RESPONSES
    CACHE_FILEPATH_UNIPROT = ""  # filepath to the file containing uniprot api queries and their final results (after processing of the url responses)
    CACHE_FILEPATH_ENSEMBL = ""  # filepath to the file containing ensembl api queries and their final results (after processing of the url responses)
    CACHE_FILEPATH_GENEONTOLOGY = ""  # filepath to the file containing gene ontology api queries and their final results (after processing of the url responses)
    CACHE_FILEPATH_GPROFILER = ""
    cached_urls = {}
    cached_uniprot = {}
    cached_ensembl = {}
    cached_geneontology = {}
    cached_gprofiler = {}
    # if any new entries have been cached - used to prevent re-saving existing caches without any changes during runtime
    delta_cached_urls = False
    delta_cached_uniprot = False
    delta_cached_ensembl = False
    delta_cached_geneontology = False
    delta_cached_gprofiler = False

    @classmethod
    def init(
        cls,
        cache_dir: str = "cache",
        store_data_atexit: bool = True,
        files_to_init:list = None
    ):
        """
        Initialises ConnectionCacher. This function must be called at the program startup in order to read
        old urls into the cls.cached_urls dictionary.

        Parameters:
          - (str) cache_dir: the cache folder
          - (bool) store_data_atexit: if True, will only store data at program exit. If False, will store data each time store_data is called.
          - (bool) use_cacher: if you wish to enable caching (setting this to 'False') will disable all Cacher functionalities
          - (list) files_to_init: The exact files which to init (files not specified will not be initialised). Possible values are: 
                                  "url", "uniprot", "ensembl", "gprofiler", "go"
        Usage:
            model = ReverseLookup.load_model("diabetes_angio_4/model_async_test.json") # make sure that model products are already computed
            Cacher.init()
            fetch_ortholog_products(refetch=True, run_async=False)

        NOTE: WARNING !! In order for the atexit storage to work, you mustn't run the Python program in VSCode in Debug mode. Run
        it in normal mode and finish the program execution with CTRL + C to test the functionality.
        """

        def faulty_cache_file_resolve(cache_filepath:str):
            # TODO:implement
            return 0
        
        should_init_urls = True
        should_init_uniprot = True
        should_init_ensembl = True
        should_init_go = True
        should_init_gprofiler = True
        files_to_init_all = ["url", "uniprot", "ensembl", "go", "gprofiler"]
        if files_to_init is not None:
            if not isinstance(files_to_init, list):
                files_to_init = [files_to_init]
                # determine which files not to initialise. Keep the files from 'files_to_init' set to True, whereas files not found in 'files_to_init' should be False
                for fa in files_to_init_all:
                    if fa not in files_to_init:
                        match fa:
                            case "url":
                                should_init_urls = False
                            case "uniprot":
                                should_init_uniprot = False
                            case "ensembl":
                                should_init_ensembl = False
                            case "go":
                                should_init_go = False
                            case "gprofiler":
                                should_init_gprofiler = False

        FileUtil.check_path(cache_dir, is_file=False)
        cls.store_data_atexit = store_data_atexit
        cls.CACHE_FILEPATH_URLS = os.path.join(cache_dir, "connection_cache.json").replace("\\", "/")
        cls.CACHE_FILEPATH_UNIPROT = os.path.join(cache_dir, "uniprot_cache.json").replace("\\", "/")
        cls.CACHE_FILEPATH_ENSEMBL = os.path.join(cache_dir, "ensembl_cache.json").replace("\\", "/")
        cls.CACHE_FILEPATH_GENEONTOLOGY = os.path.join(cache_dir, "geneontology_cache.json").replace("\\", "/")
        cls.CACHE_FILEPATH_GPROFILER = os.path.join(cache_dir, "gprofiler_cache.json").replace("\\", "/")
        FileUtil.check_paths(
            [
                cls.CACHE_FILEPATH_URLS,
                cls.CACHE_FILEPATH_GENEONTOLOGY,
                cls.CACHE_FILEPATH_UNIPROT,
                cls.CACHE_FILEPATH_ENSEMBL,
                cls.CACHE_FILEPATH_GPROFILER
            ]
        )
        if should_init_urls:
            cls.cached_urls = JsonUtil.load_json(cls.CACHE_FILEPATH_URLS)
        if should_init_uniprot:
            cls.cached_uniprot = JsonUtil.load_json(cls.CACHE_FILEPATH_UNIPROT)
        if should_init_ensembl:
            cls.cached_ensembl = JsonUtil.load_json(cls.CACHE_FILEPATH_ENSEMBL)
        if should_init_go:
            cls.cached_geneontology = JsonUtil.load_json(cls.CACHE_FILEPATH_GENEONTOLOGY)
        if should_init_gprofiler:
            cls.cached_gprofiler = JsonUtil.load_json(cls.CACHE_FILEPATH_GPROFILER)

        logger.info("Cacher load dictionary response counts:")
        logger.info(f"  - urls: {len(cls.cached_urls)}")
        logger.info(f"  - uniprot: {len(cls.cached_uniprot)}")
        logger.info(f"  - ensembl: {len(cls.cached_ensembl)}")
        logger.info(f"  - geneontology: {len(cls.cached_geneontology)}")

        if store_data_atexit:  # register the save_data function to be called on program exit
            logger.info("Register at exit save data for Cacher.")
            atexit.register(cls.save_data)

        cls.is_init = True  # to highlight that init was successful
        
        # if any new entries have been cached - used to prevent re-saving existing caches without any changes during runtime
        cls.delta_cached_urls = False
        cls.delta_cached_uniprot = False
        cls.delta_cached_ensembl = False
        cls.delta_cached_geneontology = False
        cls.delta_cached_gprofiler = False

    @classmethod
    def store_data(
        cls, data_location: str, data_key: str, data_value, timestamp: str = ""
    ):
        """
        Stores
            {
            "data_key":
                "data_value": data_value,
                "timestamp": timestamp
            }

        inside a particular json cache file, based on 'data_location'.
        The data_location options are:
          - "url": -> filepath = cache/connection_cache.json
          - "uniprot" -> filepath = cache/uniprot_cache.json
          - "ensembl" -> filepath = cache/ensembl_cache.json
          - "go" -> filepath = cache/geneontology_cache.json
          - "gprofiler" -> filepath = cache/gprofiler_cache.json

        Params:
          - (str) data_location: either 'url', 'uniprot', 'ensembl' or 'go'
          - (str) data_key: the key under which the data will be stored. The keys for function return values should be stored
                            in the format [class_name][function_name][custom_function_parameter_values]
          - (any) data_value: the value of the data being stored at data_key. Can be a whole JSON https response, or a list of processed values after function parsing of the json response etc.
          - (str) timestamp: optional, timestamps are automatically calculated inside this function if not provided

        Url storage is intended for intermediate url query responses. Consider the following url query: f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none":
        Without request caching, this is the code algorithm:

            url = ...
            response = (Session).get(url)
            response_json = response.json()

        With request caching, the code algorithm is slightly modified:

            url = ...
            previous_response = Cacher.get_data(data_location="urls", data_key=url)
            if previous_response != None:
                response_json = previous_response
            else:
                response = (Session).get(url)
                response_json = response.json()
                Cacher.store_data(data_location="urls", data_key=url, data_value=response_json)
            # process response_json

        Alternatively, in the case of storage the response jsons of the queried urls, you can use ConnectionCacher:

            url = ...
            previous_response = ConnectionCacher.get_url_response(url)
            if previous_response != None:
                response_json = previous_response
            else:
                response = (Session).get(url)
                response_json = response.json()
                ConnectionCacher.store_url(url, response=response_json)
            # process response_json

        With this code, if the algorithm encounters and already queried url, it will pull its old response,
        rather than query a new one.
        """ 
        if not hasattr(cls, "is_init"):
            logger.warning(f"MAJOR WARNING: Cacher is not storing data. Was it initialised correctly?")
            return

        if cls.is_init is False:
            cls.init()  # attempt cacher init, if the user forgot to initialise it

        cached_data = {}
        # determine type of cached data
        match data_location:
            case "url":
                cached_data = cls.cached_urls
                cls.delta_cached_urls = True
            case "uniprot":
                cached_data = cls.cached_uniprot
                cls.delta_cached_uniprot = True
            case "ensembl":
                cached_data = cls.cached_ensembl
                cls.delta_cached_ensembl = True
            case "go":
                cached_data = cls.cached_geneontology
                cls.delta_cached_geneontology = True
            case "gprofiler":
                cached_data = cls.cached_gprofiler
                cls.delta_cached_gprofiler = True

        # calculate current time
        if timestamp == "":
            timestamp = Timer.get_current_time()

        # bugfix: some urls return the following response: {'error': 'No valid lookup found for symbol Oxct2a'}
        # if this happens, do not store data
        if data_location == "url" and "error" in data_value:
            logger.warning(f"Error in data value, aborting cache store. data_value: {data_value}")
            return

        # update cached_data
        if data_key not in cached_data:
            cached_data[data_key] = {"data_value": data_value, "timestamp": timestamp}
        else:  # this data_key already exists in previous data
            previous_data_timestamp = cached_data[data_key]["timestamp"]
            if Timer.compare_time(previous_data_timestamp, timestamp) is True:  # will return true, if timestamp > previous_url_timestamp (timestamp is logged later in time than previous_url_timestamp)
                if data_value is not None:
                    cached_data[data_key] = {
                        "data_value": data_value,
                        "timestamp": timestamp,
                    }

        # update class values with new cached data and save
        if not cls.store_data_atexit:
            match data_location:
                case "url":
                    cls.cached_urls = cached_data
                    JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH_URLS)
                case "uniprot":
                    cls.cached_uniprot = cached_data
                    JsonUtil.save_json(cls.cached_uniprot, cls.CACHE_FILEPATH_UNIPROT)
                case "ensembl":
                    cls.cached_ensembl = cached_data
                    JsonUtil.save_json(cls.cached_ensembl, cls.CACHE_FILEPATH_ENSEMBL)
                case "go":
                    cls.cached_geneontology = cached_data
                    JsonUtil.save_json(cls.cached_geneontology, cls.CACHE_FILEPATH_GENEONTOLOGY)
                case "gprofiler":
                    cls.cached_gprofiler = cached_data
                    JsonUtil.save_json(cls.cached_gprofiler, cls.CACHE_FILEPATH_GPROFILER)

    @classmethod
    def get_data(cls, data_location: str, data_key: str, debug_log:bool=False):
        if not hasattr(cls, "is_init"):
            return None
        
        if cls.is_init is False:
            cls.init()

        cached_data = {}
        # determine type of cached data
        match data_location:
            case "url":
                cached_data = cls.cached_urls
            case "uniprot":
                cached_data = cls.cached_uniprot
            case "ensembl":
                cached_data = cls.cached_ensembl
            case "go":
                cached_data = cls.cached_geneontology
            case "gprofiler":
                cached_data = cls.cached_gprofiler

        if cached_data != {}:
            if data_key in cached_data:
                if debug_log:
                    logger.info(f"Successfully returned old data for {data_key}.")
                return_value = cached_data[data_key]["data_value"]

                # bugfix: some urls return the following response: {'error': 'No valid lookup found for symbol Oxct2a'}
                # if such a stored url response is read, return None
                if "error" in return_value:
                    return None
                # return_value doesn't contain error -> return it
                return return_value
            else:
                return None
        else:
            return None

    @classmethod
    def save_data(cls):
        """
        Saves 'cached_urls', 'cached_uniprot' and 'cached_ensembl' to their respective
        cache filepaths (CACHE_FILEPATH_URLS, CACHE_FILEPATH_UNIPROT, CACHE_FILEPATH_ENSEMBL)
        """
        logger.info("Cacher is saving data. Please, be patient.")
        if cls.delta_cached_urls:
            JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH_URLS)
        if cls.delta_cached_uniprot:
            JsonUtil.save_json(cls.cached_uniprot, cls.CACHE_FILEPATH_UNIPROT)
        if cls.delta_cached_ensembl:
            JsonUtil.save_json(cls.cached_ensembl, cls.CACHE_FILEPATH_ENSEMBL)
        if cls.delta_cached_geneontology:
            JsonUtil.save_json(cls.cached_geneontology, cls.CACHE_FILEPATH_GENEONTOLOGY)
        if cls.delta_cached_gprofiler:
            JsonUtil.save_json(cls.cached_gprofiler, cls.CACHE_FILEPATH_GPROFILER)
        logger.info("Successfully saved url, uniprot, ensembl, geneontology and gprofiler cache.")

    @classmethod
    def clear_cache(cls, cache_to_clear: str):
        """
        Clears the specified 'cache_to_clear', which must be one of the following:
          - "url"
          - "uniprot"
          - "ensembl"
          - "go"

        If cache_to_clear is set to "ALL", every cache will be cleared.
        """
        if cache_to_clear == "ALL":
            FileUtil.clear_file(cls.CACHE_FILEPATH_URLS, replacement_text="{}")
            FileUtil.clear_file(cls.CACHE_FILEPATH_UNIPROT, replacement_text="{}")
            FileUtil.clear_file(cls.CACHE_FILEPATH_ENSEMBL, replacement_text="{}")
            FileUtil.clear_file(cls.CACHE_FILEPATH_GENEONTOLOGY, replacement_text="{}")
            FileUtil.clear_file(cls.CACHE_FILEPATH_GPROFILER, replacement_text="{}")
            cls.cached_urls = {}
            cls.cached_uniprot = {}
            cls.cached_ensembl = {}
            cls.cached_geneontology = {}
            cls.cached_gprofiler = {}
            logger.info("Cleared entire cache.")
            return

        filepath_to_clear = ""
        match cache_to_clear:
            case "url":
                filepath_to_clear = cls.CACHE_FILEPATH_URLS
                cls.cached_urls = {}
            case "uniprot":
                filepath_to_clear = cls.CACHE_FILEPATH_UNIPROT
                cls.cached_uniprot = {}
            case "ensembl":
                filepath_to_clear = cls.CACHE_FILEPATH_ENSEMBL
                cls.cached_ensembl = {}
            case "go":
                filepath_to_clear = cls.CACHE_FILEPATH_GENEONTOLOGY
                cls.cached_geneontology = {}
            case "gprofiler":
                filepath_to_clear = cls.CACHE_FILEPATH_GPROFILER
                cls.cached_gprofiler = {}
        FileUtil.clear_file(
            filepath=filepath_to_clear, replacement_text="{}"
        )  # set empty json to cache file
        logger.info(f"Cache '{cache_to_clear}' is cleared.")


class ConnectionCacher(Cacher):
    """
    ConnectionCacher accesses root/cache/connection_cache.json in order to store and retrieve old
    url connections and their belonging data (url response, time of request)

    A newer implementation, which combines connection caching, as well as caching of processed uniprot or
    ensembl function results, is the Cacher class. It is advisable to use the Cacher class with the
    parameter "url" as data_location in place of ConnectionCacher.

    NOTE: We could make three implementations of Cacher -> ConnectionCacher, UniprotCacher, EnsemblCacher,
    but that would add too complex functionality, which can be reasonably implemented in a single class.
    """

    CACHE_FILEPATH = "cache/connection_cache.json"
    cached_urls = {}

    @classmethod
    def init(cls):
        """
        Initialises ConnectionCacher. This function must be called at the program startup in order to read
        old urls into the cls.cached_urls dictionary.
        """
        cls.CACHE_FILEPATH = "cache/connection_cache.json"
        cls.cached_urls = JsonUtil.load_json(cls.CACHE_FILEPATH)

    @classmethod
    def store_url(cls, url: str, response, timestamp: str = ""):
        """
        Stores the 'url' as the key, it's value is a dictionary comprised of 'response' and 'timestamp'.
        The key-value pair is stored in root/cache/connection_cache.json. If timestamp is not provided, then
        a timestamp will be calculated inside this function.

        Json outline:
        {url1 -> {"response": response1, "timestamp": timestamp1}},
        {url2 -> {"response": response2, "timestamp": timestamp2}},
        ...
        """
        # data = JsonUtil.load_json(cls.CACHE_FILEPATH)

        data = cls.cached_urls
        if timestamp == "":
            timestamp = Timer.get_current_time()

        if "url" not in data:  # url doesn't exist in previous data -> add it
            data[url] = {
                "response": response,
                "timestamp": timestamp,
            }  # add new element
            cls.cached_urls = data  # update cached urls
            JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH)  # save cached urls
        else:  # this url already exists in previous data
            # previous_url_response = data[url]["response"]
            previous_url_timestamp = data[url]["timestamp"]
            if (
                Timer.compare_time(previous_url_timestamp, timestamp) is True
            ):  # will return true, if timestamp > previous_url_timestamp (timestamp is logged later in time than previous_url_timestamp)
                if response is not None:
                    data[url] = {
                        "response": response,
                        "timestamp": timestamp,
                    }  # add new element
                    cls.cached_urls = data  # update cached urls
                    JsonUtil.save_json(
                        cls.cached_urls, cls.CACHE_FILEPATH
                    )  # save cached urls

    @classmethod
    def get_url_response(cls, url: str, debug_log:bool=False):
        """
        Obtains the response of the 'url' from previously cached urls, if the same url already exists.
        Previously cached urls and their responses are stored in root/cache/connection_cache.json.

        Returns None either if url doesn't exist or if the response of the url is stored as None.
        """
        if cls.cached_urls == {}:
            logger.warning(
                "Cached_urls variable is empty! Did you forget to call"
                " ConnectionCacher.init()?"
            )
            cls.init()

        if cls.cached_urls != {}:
            if url in cls.cached_urls:
                # TODO: implement url age option, eg. if the user selects "previous month", if the url is older than that, return None
                if debug_log:
                    logger.info(
                        f"Cached response for {url}: {cls.cached_urls[url]['response']}"
                    )
                return cls.cached_urls[url]["response"]
            else:  # url wasn't found
                return None
        else:
            return None

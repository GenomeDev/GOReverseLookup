from __future__ import annotations
import asyncio
import aiohttp
from contextlib import asynccontextmanager
from typing import List, Dict
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from .core.ModelSettings import ModelSettings
from .core.GOTerm import GOTerm
from .core.Product import Product
from .core.Metrics import Metrics, basic_mirna_score, miRDB60predictor
from .core.miRNA import miRNA
from .web_apis.GOApi import GOApi
from .web_apis.EnsemblApi import EnsemblApi
from .web_apis.UniprotApi import UniProtApi
from .parse.GOAFParser import GOAnnotationsFile
from .parse.OrthologParsers import HumanOrthologFinder
from .parse.OBOParser import OboParser
from .util.JsonUtil import JsonUtil, JsonToClass
from .util.Timer import Timer

import logging

# from logging import config
# config.fileConfig("logging_config.py")
logger = logging.getLogger(__name__)


class TargetProcess:
    def __init__(self, name: str, direction: str) -> None:
        """
        A class representing a target process. NOT USED CURRENTLY

        Args:
            name (str)
            direction (str): + or -
            goterms (set): a set of all goterms which are
        """
        self.name = name
        self.direction = direction


class ReverseLookup:
    def __init__(
        self,
        goterms: List[GOTerm],
        target_processes: List[Dict[str, str]],
        products: List[Product] = [],
        miRNAs: List[miRNA] = [],
        miRNA_overlap_treshold: float = 0.6,
        execution_times: dict = {},
        statistically_relevant_products={},
        go_categories: List[str] = [
            "biological_process",
            "molecular_activity",
            "cellular_component",
        ],
        model_settings: ModelSettings = None,
        obo_parser: OboParser = None,
        datafile_paths: dict = {},
    ):
        """
        A class representing a reverse lookup for gene products and their associated Gene Ontology terms.

        Args:
            goterms (set): A set of GOTerm objects.
            target_processes (list): A list of dictionaries containing process names and directions.
            products (set, optional): A set of Product objects. Defaults to an empty set.
            miRNAs
            miRNA_overlap_threshold
            execution_times: a dictionary of function execution times, no use to the user. Used during model saving and loadup.
            statistically_relevant_products: TODO - load and save the model after perform_statistical_analysis is computed
            go_categories: When querying GO Term data, which GO categories should be allowed. By default, all three categories are allowed ("biological_process", "molecular_activity", "cellular_component").
                        Choosing the correct categories affects primarily the Fisher scoring, when GO Terms are queried for each product either from the GOAF or from the web. Excluding some GO categories (such as cellular_component)
                        when researching only GO Terms connected to biological processes and/or molecular activities helps to produce more accurate statistical scores.
            model_settings: used for saving and loading model settings
            obo_parser: Used for parsing the Gene Ontology's .obo file. If it isn't supplied, it will be automatically created
            datafile_paths: The paths to 3rd party database files. This is a dictionary, each element holds the data file type as key and the data file path as the value.
                            Datafile paths are queried using (ReverseLookup).get_datafile_path(datafile_type) function, where allowed datafile_types are:
                            'go_obo_filepath', 'goaf_filepath', 'zfin_human_ortho_mapping_filepath', 'mgi_human_ortho_mapping_filepath', 'rgd_human_ortho_mapping_filepath', 'xenbase_human_ortho_mapping_filepath'
        """
        self.goterms = goterms
        self.products = products
        self.target_processes = target_processes
        self.miRNAs = miRNAs
        self.miRNA_overlap_treshold = miRNA_overlap_treshold
        self.model_settings = model_settings

        # to enable backwards compatibility
        if datafile_paths == {}:
            logger.info("NOTE: Model datafile paths are all None!")
            datafile_paths = {
                "go_obo_filepath": None,
                "goaf_filepath": None,
                "zfin_human_ortho_mapping_filepath": None,
                "mgi_human_ortho_mapping_filepath": None,
                "rgd_human_ortho_mapping_filepath": None,
                "xenbase_human_ortho_mapping_filepath": None,
            }
        else:
            self.datafile_paths = datafile_paths

        self.execution_times = (
            execution_times  # dict of execution times, logs of runtime for functions
        )
        self.timer = Timer()

        # placeholder to populate after perform_statistical_analysis is called
        self.statistically_relevant_products = statistically_relevant_products

        # GO categories - determines which categories of GO terms are chosen during the Fisher score computation (either from the GO Annotations File or from a gene id to GO Terms API query)
        self.go_categories = go_categories
        self.goaf = GOAnnotationsFile(
            filepath=self.datafile_paths["goaf_filepath"], go_categories=go_categories
        )
        if self.goaf is None:
            logger.warning("MODEL COULD NOT CREATE A GO ANNOTATIONS FILE!")
            logger.warning(
                f"  - goaf_filepath = {self.datafile_paths['goaf_filepath']}"
            )

        self.go_api = GOApi() # this enables us to use go_api inside Metrics.py, as importing GOApi inside Metrics.py creates circular imports.

        if obo_parser is not None:
            self.obo_parser = obo_parser
        else:
            # obo parser is non
            if (
                self.datafile_paths["go_obo_filepath"] is not None
                and self.datafile_paths["go_obo_filepath"] != ""
            ):
                self.obo_parser = OboParser(
                    obo_filepath=self.datafile_paths["go_obo_filepath"]
                )
            else:
                self.obo_parser = OboParser()

    def set_model_settings(self, model_settings: ModelSettings):
        """
        Sets self.model_settings to the model settings supplied in the parameter.
        """
        self.model_settings = model_settings
    
    def get_goterm(self, goterm_id:str):
        """
        Returns the GO Term object associated with the input 'goterm_id' string
        """
        for goterm in self.goterms:
            if goterm.id == goterm_id:
                return goterm

    def set_model_setting(self, setting: str, value):
        """
        If the attribute 'setting' in self.model_settings (ModelSettings) exists, sets its value
        to 'value'
        """
        if hasattr(self.model_settings, setting):
            setattr(self.model_settings, setting, value)
        else:
            logger.warning(f"ModelSettings object has no attribute {setting}!")

    def get_datafile_path(self, datafile_type: str):
        """
        Gets the path to a data file specified by datafile_type. Datafile paths are loaded from the 'filepaths' section in input.txt.
        Allowed datafile types (the values of datafile_type parameter) are:
          - 'go_obo_filepath'
          - 'goaf_filepath'
          - 'zfin_human_ortho_mapping_filepath'
          - 'mgi_human_ortho_mapping_filepath'
          - 'rgd_human_ortho_mapping_filepath'
          - 'xenbase_human_ortho_mapping_filepath'
        """
        if self.datafile_paths != []:
            if datafile_type in self.datafile_paths:
                return self.datafile_paths[datafile_type]
            else:
                logger.warning(
                    f"Datafile type {datafile_type} doesn't exist in"
                    " self.datafile_paths. Returning None."
                )
                return None
        logger.warning(
            "self.datafile_paths wasn't initialised (is empty). Returning None."
        )
        return None

    def get_datafile_paths(self, *datafile_types: str):
        """
        Gets the paths to multiple datafile_types. Datafile paths are loaded from the 'filepaths' section in input.txt.
        Allowed datafile types (the values of datafile_type parameter) are:
          - 'go_obo_filepath'
          - 'goaf_filepath'
          - 'zfin_human_ortho_mapping_filepath'
          - 'mgi_human_ortho_mapping_filepath'
          - 'rgd_human_ortho_mapping_filepath'
          - 'xenbase_human_ortho_mapping_filepath'

        Returns a dictionary with the keys corresponding to input datafile_types and the values to the found filepaths.

        Example call:
        (ReverseLookup).get_datafile_paths('go_obo_filepath', 'goaf_filepath')
        -> returns:
            {
                'go_obo_filepath': "PATH_TO_GO_OBO_FILE" or None,
                'goaf_filepath': "PATH_TO_GOAF_FILE" or None
            }

        Note: You can call this function with only one string parameter "ALL" in order to receive a dictionary of all possible
        datafile types.
        """
        all_keys = [
            "go_obo_filepath",
            "goaf_filepath",
            "zfin_human_ortho_mapping_filepath",
            "mgi_human_ortho_mapping_filepath",
            "rgd_human_ortho_mapping_filepath",
            "xenbase_human_ortho_mapping_filepath",
        ]
        if len(datafile_types) == 1 and datafile_types[0] == "ALL":
            input_keys = all_keys
        else:
            input_keys = datafile_types

        result = {}
        for datafile_type in input_keys:
            result[datafile_type] = self.get_datafile_path(datafile_type)
        return result

    def fetch_all_go_term_names_descriptions(self, run_async=True, req_delay=0.1, max_connections=50):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_name_description method for each object.
        """
        self.timer.set_start_time()
        api = GOApi()

        if run_async is True:
            asyncio.run(
                self._fetch_all_go_term_names_descriptions_async(
                    api, req_delay=req_delay, max_connections=max_connections
                )
            )
        else:
            logger.info("Fetching GO term names and their descriptions.")
            # TODO: tqdm prevents any logger.info to be printed to console
            # tqdm.write(f"Fetching GO term names and their descriptions.")
            with logging_redirect_tqdm():
                for goterm in tqdm(self.goterms, desc="Fetch term names and descs"):
                    if (
                        goterm.name is None or goterm.description is None
                    ):  # if goterm.name or description don't exist, then attempt fetch
                        goterm.fetch_name_description(api)

        if (
            "fetch_all_go_term_names_descriptions" not in self.execution_times
        ):  # to prevent overwriting on additional runs of the same model name
            self.execution_times[
                "fetch_all_go_term_names_descriptions"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_all_go_term_names_descriptions_async(
        self, api: GOApi, req_delay=0.1, max_connections=50
    ):
        """
        Call fetch_all_go_term_names_descriptions with run_async == True to run this code.
        """
        connector = aiohttp.TCPConnector(
            limit=max_connections, limit_per_host=max_connections
        ) 
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = []
            for goterm in self.goterms:
                if goterm.name is None or goterm.description is None:
                    task = asyncio.create_task(
                        goterm.fetch_name_description_async(api=api, session=session, req_delay=req_delay)
                    )
                    tasks.append(task)
            await asyncio.gather(*tasks)

    def fetch_all_go_term_products(
        self,
        web_download: bool = True,
        run_async=True,
        recalculate: bool = False,
        delay: float = 0.2,
        run_async_options: str = "v3",
        request_params={"rows": 50000, "taxon": "NCBITaxon:9606"},
        max_connections=60,
    ):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_products method for each object.

        Args:
          - (bool) recalculate: if set to True, will recalculate (fetch again) the term's products even if they already exist (perhaps from a model loaded from data.json)
          - (bool) run_async: if True, will run web downloads asynchronously, which is 1000x faster than synchronous web downloads
          - (bool) web_download: if set to True, then the products will be downloaded using https api queries. If set to False, then the products for GO Terms will be
                                 parsed from a GOAnnotationFile (http://current.geneontology.org/products/pages/downloads.html).
          - (float) delay: the delay between async requests
          - (str) run_async_options: either v1 or v2 (for development purposes)
                - v1 created as many ClientSession objects as there are goterms -> there is no control
                  over the amount of requests sent to the server, since each ClientSession is sending only one request, but they simultaneously clutter the server.
                  The result are 504 bad gateway requests
                - v2 creates only one ClientSession object for all goterms (further divisions could be possible for maybe 2 or 4 ClientSessions to segment the goterms),
                  which allows us to control the amount of requests sent to the server. The result is that the server doesn't detect us as bots and doesn't block our requests.
                  v2 should be used.
                - v3 is the best working function and should be always used.

        Developer explanation for v1, v2 and v3 versions of async:
          - *** async version 1 ***
            code:
                tasks = []
                api = GOApi()
                for goterm in self.goterms:
                    if goterm.products == [] or recalculate == True:
                        task = asyncio.create_task(goterm.fetch_products_async_v1(api, delay=delay))
                            --- ---
                            await asyncio.sleep(delay)
                            # products = await api.get_products_async(self.id)
                            products = await api.get_products_async_notimeout(self.id)
                                --- ---
                                url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
                                connector = aiohttp.TCPConnector(limit=20, limit_per_host=20)
                                async with aiohttp.ClientSession(connector=connector) as session:
                                    response = await session.get(url, params=params)
                                    ...
                                --- ---
                            --- ---
                        tasks.append(task)
                await asyncio.gather(*tasks)

            explanation:
                The connector object is created for each GO Term object. There is no single "master" connector,
                hence connections to the server aren't controlled. The server is overloaded with connections and blocks incoming
                connections.

          - *** async version 2 ***
            code:
                api = GOApi()
                connector = aiohttp.TCPConnector(limit=max_connections,limit_per_host=max_connections) # default is 100
                async with aiohttp.ClientSession(connector=connector) as session:
                    for goterm in self.goterms:
                        url = api.get_products(goterm.id,get_url_only=True, request_params=request_params)
                        await asyncio.sleep(req_delay) # request delay
                        response = await session.get(url)
                        ...

            explanation:
                In contrast to v1, this code uses a master connector for the ClientSession, but is much slower, as the
                requests are actually sent synchronously (each response is awaited inside the for loop). Thus, this function
                doesn't achieve the purpose of async requests, but demonstrates how to limit server connections using a master ClientSession.

          - *** async version 3 ***
            code:
                connector = aiohttp.TCPConnector(limit=max_connections,limit_per_host=max_connections) # default is 100
                async with aiohttp.ClientSession(connector=connector) as session:
                    tasks = []
                    for goterm in self.goterms:
                        task = goterm.fetch_products_async_v3(session, request_params=request_params, req_delay=req_delay)
                            --- ---
                            url = f"http://api.geneontology.org/api/bioentity/function/{self.id}/genes"
                            asyncio.sleep(req_delay)
                            response = await session.get(url, params=params)
                            ...
                            --- ---
                        tasks.append(task)
                    # perform multiple tasks at once asynchronously
                    await asyncio.gather(*tasks)

            explanation:
                The v3 version of the code uses asyncio.gather, which concurrently runs the list of awaitable
                objects in the supplied parameter list. First, all execution tasks are gathered in a list, which is
                then supplied to asyncio.gather. The code also uses a master ClientSession with a custom TCPConnector object,
                which limits the maximum server connections.
        """
        logger.info("Started fetching all GO Term products.")
        self.timer.set_start_time()

        if web_download is True:
            source = GOApi()
        else:
            source = GOAnnotationsFile()

        if run_async is True:
            if run_async_options == "v1":
                asyncio.run(
                    self._fetch_all_go_term_products_async_v1(
                        recalculate=False, delay=delay
                    )
                )
            elif run_async_options == "v2":
                asyncio.run(
                    self._fetch_all_goterm_products_async_v2(
                        max_connections=max_connections,
                        request_params=request_params,
                        req_delay=delay,
                    )
                )
            elif run_async_options == "v3":
                asyncio.run(
                    self._fetch_all_goterm_products_async_v3(
                        max_connections=max_connections,
                        request_params=request_params,
                        req_delay=delay,
                    )
                )
        else:
            with logging_redirect_tqdm():
                for goterm in tqdm(self.goterms, desc="Fetch term products"):
                    if (
                        goterm.products == [] or recalculate is True
                    ):  # to prevent recalculation of products if they are already computed
                        goterm.fetch_products(source)

        # api = GOApi()
        # with logging_redirect_tqdm():
        #     for goterm in tqdm(self.goterms):
        #         goterm.fetch_products(api)

        if "fetch_all_go_term_products" not in self.execution_times:
            self.execution_times[
                "fetch_all_go_term_products"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_all_go_term_products_async_v1(
        self, recalculate: bool = False, delay: float = 0.0
    ):
        """
        Asynchronously queries the products for all GOTerm objects. Must be a web download.
        This function is 1000x times faster than it's synchronous 'fetch_all_go_term_products' counterpart

        To call this function, call 'fetch_all_go_term_products' with run_async = True [TODO]

        Args:
          - (bool) recalculate: if set to True, will recalculate (fetch again) the term's products even if they already exist (perhaps from a model loaded from data.json)
          - (float) delay: the delay between asyncio requests
        """
        tasks = []
        api = GOApi()
        for goterm in self.goterms:
            if goterm.products == [] or recalculate is True:
                # sleeping here doesnt fix the server blocking issue!
                task = asyncio.create_task(
                    goterm.fetch_products_async_v1(api, delay=delay)
                )
                tasks.append(task)
        await asyncio.gather(*tasks)

    async def _fetch_all_goterm_products_async_v2(
        self, max_connections=100, request_params={"rows": 20000}, req_delay=0.5
    ):
        """
        In comparison to _fetch_all_go_term_products_async, this function doesn't overload the server and cause the server to block our requests.
        """
        APPROVED_DATABASES = [
            ["UniProtKB", ["NCBITaxon:9606"]],
            ["ZFIN", ["NCBITaxon:7955"]],
            # ["RNAcentral", ["NCBITaxon:9606"]],
            ["Xenbase", ["NCBITaxon:8364"]],
            ["MGI", ["NCBITaxon:10090"]],
            ["RGD", ["NCBITaxon:10116"]],
        ]
        api = GOApi()

        connector = aiohttp.TCPConnector(
            limit=max_connections, limit_per_host=max_connections
        )  # default is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            for goterm in self.goterms:
                url = api.get_products(
                    goterm.id, get_url_only=True, request_params=request_params
                )
                await asyncio.sleep(req_delay)  # request delay
                response = await session.get(url)

                if (
                    response.status != 200
                ):  # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                    logger.warning(
                        f"HTTP Error when parsing {goterm.id}. Response status ="
                        f" {response.status}"
                    )
                    goterm.http_error_codes["products"] = (
                        f"HTTP Error: status = {response.status}, reason ="
                        f" {response.reason}"
                    )

                data = await response.json()
                products_set = set()
                for assoc in data["associations"]:
                    if assoc["object"]["id"] == goterm.id and any(
                        (
                            database[0] in assoc["subject"]["id"]
                            and any(
                                taxon in assoc["subject"]["taxon"]["id"]
                                for taxon in database[1]
                            )
                        )
                        for database in APPROVED_DATABASES
                    ):
                        product_id = assoc["subject"]["id"]
                        products_set.add(product_id)

                products = list(products_set)
                logger.info(f"Fetched products for GO term {goterm.id}")
                goterm.products = products

    # IMPROVE SPEED UP USING ASYNCIO.GATHER: Instead of awaiting each request individually in a loop, you can use asyncio.gather()
    # to concurrently execute multiple requests. This allows the requests to be made in parallel, which can significantly improve performance.
    async def _fetch_all_goterm_products_async_v3(
        self,
        max_connections=100,
        request_params={"rows": 50000, "taxon": "NCBITaxon:9606"},
        req_delay=0.5,
        recalculate: bool = False,
    ):
        """
        In comparison to (GOApi)._fetch_all_go_term_products_async, this function doesn't overload the server and cause the server to block our requests.
        In comparison to the v2 version of this function (inside GOApi), v3 uses asyncio.gather, which speeds up the async requests.
        """
        connector = aiohttp.TCPConnector(
            limit=max_connections, limit_per_host=max_connections
        )  # default is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = []
            for goterm in self.goterms:
                if goterm.products == [] or recalculate is True:
                    task = goterm.fetch_products_async_v3(
                        session, request_params=request_params, req_delay=req_delay
                    )
                    tasks.append(task)
            # perform multiple tasks at once asynchronously
            await asyncio.gather(*tasks)

    def create_products_from_goterms(self) -> None:
        """
        This method creates Product objects from the set of products contained in each GOTerm object and
        adds them to the ReverseLookup object's products list.

        The method iterates over each GOTerm object in the goterms set and retrieves the set of products associated
        with that GOTerm object. It then adds these products to a products_set, which is a set object that ensures
        that no duplicate products are added.

        Finally, the method iterates over each product in the products_set and creates a new Product object from the
        product ID using the Product.from_dict() classmethod. The resulting Product objects are added to the
        ReverseLookup object's products list.

        Args:
            None

        Returns:
            None
        """
        logger.info(
            f"Creating products from GO Terms. Num goterms = {len(self.goterms)}"
        )
        self.timer.set_start_time()

        # Create an empty set to store unique products
        products_set = set()

        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            products_set.update(term.products)

        # Iterate over each product in the products_set and create a new Product object from the product ID using the
        # Product.from_dict() classmethod. Add the resulting Product objects to the ReverseLookup object's products list.
        for product in products_set:
            if ":" in product:
                self.products.append(Product.from_dict({"id_synonyms": [product]}))
            else:
                self.products.append(
                    Product.from_dict({"id_synonyms": [product], "genename": product})
                )
        logger.info(f"Created {len(self.products)} Product objects from GOTerm object definitions")

        if "create_products_from_goterms" not in self.execution_times:
            self.execution_times[
                "create_products_from_goterms"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

        def check_exists(product_id: str) -> bool:
            """
            Checks if the product_id already exists among self.products. When loading using ReverseLookup.load_model(data.json),
            check_exists has to be used in order to prevent product duplications.

            Returns: True, if product_id already exists in self.products
            """
            for existing_product in self.products:
                if product_id in existing_product.id_synonyms:
                    return True
            return False

        logger.info("Creating products from GO Terms")
        self.timer.set_start_time()

        # Create an empty set to store unique products
        products_set = set()

        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            products_set.update(term.products)

        # Iterate over each product in the products_set and create a new Product object from the product ID using the
        # Product.from_dict() classmethod. Add the resulting Product objects to the ReverseLookup object's products list.
        i = 0
        for (
            product
        ) in products_set:  # here, each product is a product id, eg. 'MGI:1343124'
            if check_exists(product) is True:
                continue
            if ":" in product:
                self.products.append(Product.from_dict({"id_synonyms": [product]}))
            else:
                self.products.append(
                    Product.from_dict({"id_synonyms": [product], "genename": product})
                )
            i += 1

        logger.info(f"Created {i} Product objects from GOTerm object definitions")

        if "create_products_from_goterms" not in self.execution_times:
            self.execution_times[
                "create_products_from_goterms"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def fetch_ortholog_products(
        self,
        refetch: bool = False,
        run_async=True,
        use_goaf=False,
        max_connections=15,
        req_delay=0.1,
        semaphore_connections=5,
    ) -> None:
        """
        This function tries to find the orthologs to any non-uniprot genes (products) associated with a GO Term.

        Args:
          - (bool) refetch: if True, will fetch the ortholog products for all Product instances again, even if some Product instances already have their orthologs fetched.
          - (bool) run_async: if True, will send requests asynchronously
          - (int) max_connections: the maximum amount of connections the asynchronous client session will send to the server
          - (float) req_delay: the delay between connections in secondsž

        This function relies on request caching. It will cache the http requests to the server. When using async requests with ortholog fetch, the first full run of all products is successful, but if
        the user decides to run async ortholog query for the same products again, the server will start sending 429:TooManyRequests error. Therefore, this option saves (caches)
        the requests in a dictionary, where the key is the request url and the value is a dictionary with request info (response, query timestamp, ...). When using async requests,
        the responses of the previous cached requests are used, if the request urls are the same. TODO: daj userju možnost, da selecta "starost" requesta aka da lahko v funkcijo poslje "7 dni"
        in bo potem uporabilo requeste, ki so "mlajsi" od 7 dni.

        NOTE: This function is recalculation-optimised based on the "genename" field of the Product. If the model is loaded from data.json and a specific
        Product already had orthologs fetched, then it is skipped during the fetch_ortholog call.

        When fetching products (genes / gene products) from Gene Ontology for a specific GO Term:
            (GOTerm).fetch_products()

            api = GOApi()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_products(api)

        The resulting products can be from any of the following databases: UniProtKB, ZFIN, Xenbase, MGI, RGD. For subsequent
        Homo-Sapiens-only product analysis, it is important to find, if human ortholog genes exist for the products, fetched from a non-uniprot
        databases.

        Usage and calling:
            products = ... # define a list of Product instances
            human_ortholog_finder = HumanOrthologFinder()
            uniprot_api = UniProtAPI()
            ensembl_api = EnsemblAPI()

            for product in products:
                product.fetch_ortholog(human_ortholog_finder, uniprot_api, ensembl_api)

        """
        logger.info("Started fetching ortholog products.")
        self.timer.set_start_time()

        try:
            if use_goaf is True:
                """
                Use the GO Annotations File to query orthologs.
                # TODO: implement async goaf parsing
                """
                # TODO: with logging_redirect_tqdm
                # TODO: remove this and use_goaf

            elif run_async is True:
                asyncio.run(
                    self._fetch_ortholog_products_async(
                        refetch=refetch,
                        max_connections=max_connections,
                        req_delay=req_delay,
                        semaphore_connections=semaphore_connections,
                    )
                )
            else:
                third_party_db_files = self.get_datafile_paths("ALL")
                human_ortholog_finder = HumanOrthologFinder(
                    goaf=self.goaf,
                    zfin_filepath=third_party_db_files[
                        "zfin_human_ortho_mapping_filepath"
                    ],
                    xenbase_filepath=third_party_db_files[
                        "xenbase_human_ortho_mapping_filepath"
                    ],
                    mgi_filepath=third_party_db_files[
                        "mgi_human_ortho_mapping_filepath"
                    ],
                    rgd_filepath=third_party_db_files[
                        "rgd_human_ortho_mapping_filepath"
                    ],
                )
                uniprot_api = UniProtApi()
                ensembl_api = EnsemblApi()

                with logging_redirect_tqdm():
                    for product in tqdm(
                        self.products, desc="Fetch ortholog products"
                    ):  # Iterate over each Product object in the ReverseLookup object.
                        # Check if the Product object doesn't have a UniProt ID or genename or ensg_id -> these indicate no ortholog computation has been performed yet
                        # if product.genename == None or refetch == True: # product.genename was still None for a lot of products, despite calling fetch_orthologs
                        if product.had_orthologs_computed is False or refetch is True:
                            # If it doesn't, fetch UniProt data for the Product object.
                            product.fetch_ortholog(
                                human_ortholog_finder,
                                uniprot_api,
                                ensembl_api,
                                goaf=self.goaf,
                            )
                            product.had_orthologs_computed = True
        except Exception as e:
            # If there was an exception while fetching UniProt data, save all the Product objects to a JSON file.
            self.save_model("crash_products.json")
            # Re-raise the exception so that the caller of the method can handle it.
            raise e

        if "fetch_ortholog_products" not in self.execution_times:
            self.execution_times[
                "fetch_ortholog_products"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_ortholog_products_async(
        self,
        refetch: bool = True,
        max_connections=100,
        req_delay=0.5,
        semaphore_connections=10,
    ):
        """
        code: [TODO: delete this after testing is done]
        connector = aiohttp.TCPConnector(limit=max_connections,limit_per_host=max_connections) # default is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = []
            for goterm in self.goterms:
                task = goterm.fetch_products_async_v3(session, request_params=request_params, req_delay=req_delay)
                    --- ---
                    url = f"http://api.geneontology.org/api/bioentity/function/{self.id}/genes"
                    asyncio.sleep(req_delay)
                    response = await session.get(url, params=params)
                    ...
                    --- ---
                tasks.append(task)
            # perform multiple tasks at once asynchronously
            await asyncio.gather(*tasks)
        """

        @asynccontextmanager
        async def create_session():
            connector = aiohttp.TCPConnector(
                limit=max_connections, limit_per_host=max_connections
            )
            session = aiohttp.ClientSession(connector=connector)
            try:
                yield session
            finally:
                await session.close()

        third_party_db_files = self.get_datafile_paths("ALL")
        human_ortholog_finder = HumanOrthologFinder(
            goaf=self.goaf,
            zfin_filepath=third_party_db_files["zfin_human_ortho_mapping_filepath"],
            xenbase_filepath=third_party_db_files["xenbase_human_ortho_mapping_filepath"],
            mgi_filepath=third_party_db_files["mgi_human_ortho_mapping_filepath"],
            rgd_filepath=third_party_db_files["rgd_human_ortho_mapping_filepath"],
        )
        uniprot_api = UniProtApi()
        ensembl_api = EnsemblApi()
        ensembl_api.async_request_sleep_delay = req_delay
        uniprot_api.async_request_sleep_delay = req_delay
        # goaf = GOAnnotationsFile(third_party_db_files['goaf_filepath'])

        connector = aiohttp.TCPConnector(
            limit=max_connections, limit_per_host=max_connections
        )
        semaphore = asyncio.Semaphore(semaphore_connections)
        async with aiohttp.ClientSession(connector=connector) as session:
        # async with create_session() as session:
            tasks = []
            for product in self.products:
                if product.had_orthologs_computed is False or refetch is True:
                    # task = product.fetch_ortholog_async(session, human_ortholog_finder, uniprot_api, ensembl_api)
                    task = product.fetch_ortholog_async_semaphore(
                        session,
                        semaphore,
                        self.goaf,
                        human_ortholog_finder,
                        uniprot_api,
                        ensembl_api,
                    )
                    tasks.append(task)
                    product.had_orthologs_computed = True
            await asyncio.gather(*tasks)

        logger.info(
            "During ortholog query, there were"
            f" {len(ensembl_api.ortholog_query_exceptions)} ensembl api exceptions and"
            f" {len(uniprot_api.uniprot_query_exceptions)} uniprot api exceptions."
        )

        # logger.debug(f"Printing exceptions:")
        # i = 0
        # for exception_dict in ensembl_api.ortholog_query_exceptions:
        #    product_id = exception_dict.keys()[0]
        #    exception = exception_dict[product_id]
        #    logger.debug(f"[{i}] :: {product_id} : {exception}")
        #    i += 1

    def prune_products(self) -> None:
        logger.info("Started pruning products.")
        self.timer.set_start_time()

        # Create a dictionary that maps genename to a list of products
        reverse_genename_products = {}
        for product in self.products:
            if product.genename is not None:
                reverse_genename_products.setdefault(product.genename, []).append(
                    product
                )

        # For each ENSG that has more than one product associated with it, create a new product with all the synonyms
        # and remove the individual products from the list
        for genename, product_list in reverse_genename_products.items():
            if len(product_list) > 1:
                id_synonyms = []
                for product in product_list:
                    self.products.remove(product)
                    id_synonyms.extend(product.id_synonyms)
                # Create a new product with the collected information and add it to the product list
                self.products.append(
                    Product(
                        id_synonyms,
                        product_list[0].genename,
                        product_list[0].uniprot_id,
                        product_list[0].description,
                        product_list[0].ensg_id,
                        product_list[0].enst_id,
                        product_list[0].refseq_nt_id,
                        product_list[0].mRNA,
                        {},
                        product_list[0].had_orthologs_computed,
                        product_list[0].had_fetch_info_computed,
                    )
                )

        if "prune_products" not in self.execution_times:
            self.execution_times["prune_products"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def fetch_product_infos(
        self,
        refetch: bool = False,
        run_async=True,
        max_connections=15,
        semaphore_connections=5,
        req_delay=0.1,
        required_keys=["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"],
    ) -> None:
        # TODO: ensembl support batch request

        logger.info("Started fetching product infos.")
        self.timer.set_start_time()

        if run_async:
            # async mode
            asyncio.run(
                self._fetch_product_infos_async(
                    required_keys=required_keys,
                    refetch=refetch,
                    max_connections=max_connections,
                    req_delay=req_delay,
                    semaphore_connections=semaphore_connections,
                )
            )
        else:
            # sync mode
            uniprot_api = UniProtApi()
            ensembl_api = EnsemblApi()
            try:
                # Iterate over each Product object in the ReverseLookup object.
                with logging_redirect_tqdm():
                    for product in tqdm(self.products, desc="Fetch product infos"):
                        # Check if the Product object doesn't have a UniProt ID.
                        # if any(attr is None for attr in [product.genename, product.description, product.enst_id, product.ensg_id, product.refseq_nt_id]) and (product.uniprot_id or product.genename or product.ensg_id): # some were still uninitialised, despite calling fetch_product_infos
                        if product.had_fetch_info_computed is False or refetch is True:
                            # If it doesn't, fetch UniProt data for the Product object.
                            product.fetch_info(
                                uniprot_api, ensembl_api, required_keys=required_keys
                            )
                            product.had_fetch_info_computed = True
                            if product.had_fetch_info_computed is False:
                                logger.warning(
                                    "had_fetch_info_computed IS FALSE despite being"
                                    f" called for {product.id_synonyms}, genename ="
                                    f" {product.genename}"
                                )
            except Exception as e:
                raise e

        if "fetch_product_infos" not in self.execution_times:
            self.execution_times[
                "fetch_product_infos"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_product_infos_async(
        self,
        required_keys=["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"],
        refetch: bool = False,
        max_connections=50,
        req_delay=0.1,
        semaphore_connections=5,
    ):
        uniprot_api = UniProtApi()
        ensembl_api = EnsemblApi()
        uniprot_api.async_request_sleep_delay = req_delay
        ensembl_api.async_request_sleep_delay = req_delay

        connector = aiohttp.TCPConnector(
            limit=max_connections, limit_per_host=max_connections
        )
        semaphore = asyncio.Semaphore(semaphore_connections)

        async with aiohttp.ClientSession(connector=connector) as session:
            # async with create_session() as session:
            tasks = []
            for product in self.products:
                if product.had_fetch_info_computed is False or refetch is True:
                    # task = product.fetch_ortholog_async(session, human_ortholog_finder, uniprot_api, ensembl_api)
                    task = product.fetch_info_async_semaphore(
                        session, semaphore, uniprot_api, ensembl_api, required_keys
                    )
                    tasks.append(task)
                    product.had_fetch_info_computed = True
            await asyncio.gather(*tasks)

    def score_products(
        self, score_classes: List[Metrics], recalculate: bool = True
    ) -> None:
        """
        Scores the products of the current ReverseLookup model. This function allows you to pass a custom or a pre-defined scoring algorithm,
        which is of 'Metrics' type (look in Metrics.py), or a list of scoring algorithms. Each Product class of the current ReverseLookup instance products (self.products)
        has a member field 'scores'. For each product, score is computed and saved to the product's 'scores' dictionary as a mapping between the
        scoring algorithm's name (eg. "adv_score") and the corresponding product's score computed with this scoring algorithm (eg. 14.6).
        If multiple scoring algorithms are used, then the product's 'scores' dictionary will have multiple elements, each a mapping between
        the scoring algorithm's name and the corresponding score.

        Note: if a miRNA scoring algorithm is passed, such as 'basic_miRNA_score', this function redirects to self.score_miRNAs(...)

        Parameters:
          - score_classes: A subclass (implementation) of the Metrics superclass (interface). Current pre-defined Metrics implementations subclasses
                         are 'adv_product_score', 'nterms', 'inhibited_products_id', 'basic_mirna_score'.
          - (bool) recalculate: if True, will recalculate scores if they already exist. If False, will skip recalculations.

        Calling example:
        (1) Construct a ReverseLookup model
        model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

        (2) Create one or more Metrics scoring implementations for the model:
        adv_score = adv_product_score(model)
        nterms_score = nterms(model)

        (3) Call the score_products on the model using the Metrics scoring implementations
        model.score_products([adv_score, nterms_score])
        """
        logger.info("Started scoring products.")
        self.timer.set_start_time()

        if not isinstance(score_classes, list):
            score_classes = [score_classes]

        # perform scoring of each product (gene)
        with logging_redirect_tqdm():
            # iterate over each Product object in self.products and score them using the Scoring object
            for product in tqdm(
                self.products, "Scoring products"
            ):  # each Product has a field scores - a dictionary between a name of the scoring algorithm and it's corresponding score
                for _score_class in score_classes:
                    # NOTE: Current miRNA scoring (self.score_miRNAs) performs miRNA scoring holistically - in one call for all miRNAs in self.miRNAs. It is pointless to call this function here, as it needs to
                    # be called only once. Here, a function for miRNA scoring has to be called, which displays the top N miRNAs, which bind to the specific product.
                    #
                    # if isinstance(_score_class, basic_mirna_score):
                    #    self.score_miRNAs(_score_class, recalculate=recalculate)
                    #    continue
                    if isinstance(_score_class, basic_mirna_score):
                        # just continue, see explanation above
                        continue

                    if (
                        _score_class.name in product.scores and recalculate is True
                    ):  # if score already exists and recalculate is set to True
                        product.scores[_score_class.name] = _score_class.metric(
                            product
                        )  # create a dictionary between the scoring algorithm name and it's score for current product
                    elif (
                        _score_class.name not in product.scores
                    ):  # if score doesn't exist yet
                        product.scores[_score_class.name] = _score_class.metric(product)

        # calculate Benjamini-Hochberg FDR correction
        for _score_class in score_classes:
            if isinstance(_score_class, basic_mirna_score):
                # score miRNAs holistically here, see # NOTE
                self.score_miRNAs(_score_class, recalculate=recalculate)
                continue

            i = 0
            p_values = []
            if (
                _score_class.name == "fisher_test"
                or _score_class.name == "binomial_test"
            ):
                for product in self.products:
                    for process in self.target_processes:
                        for direction in ["+", "-"]:
                            if (
                                "error"
                                in product.scores[_score_class.name][
                                    f"{process['process']}{direction}"
                                ]
                            ):  # check if there is "error" key
                                continue
                            p_values.append(
                                product.scores[_score_class.name][
                                    f"{process['process']}{direction}"
                                ]["pvalue"]
                            )
                # apply Benjamini-Hochberg FDR correction
                if len(p_values) > 0:
                    from statsmodels.stats.multitest import multipletests

                    reject, p_corrected, _, _ = multipletests(
                        p_values, alpha=0.05, method="fdr_bh"
                    )
                    for product in self.products:
                        for process in self.target_processes:
                            for direction in ["+", "-"]:
                                if (
                                    "error"
                                    in product.scores[_score_class.name][
                                        f"{process['process']}{direction}"
                                    ]
                                ):  # check if there is "error" key
                                    continue
                                product.scores[_score_class.name][
                                    f"{process['process']}{direction}"
                                ]["pvalue_corr"] = p_corrected[i]
                                i += 1

        if "score_products" not in self.execution_times:
            self.execution_times["score_products"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def fetch_mRNA_sequences(self, refetch=False) -> None:
        logger.info("Started fetching mRNA sequences.")
        self.timer.set_start_time()

        try:
            ensembl_api = EnsemblApi()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products, desc="Fetch mRNA seqs"):
                    # Check if the Product object doesn't have a EnsemblID
                    if (
                        product.mRNA == -1 and refetch is False
                    ):  # product mRNA was already fetched, but unsuccessfully
                        continue
                    if product.mRNA is None and product.enst_id is not None:
                        # If it has, fetch mRNA sequence data for the Product object.
                        product.fetch_mRNA_sequence(ensembl_api)
        except Exception as e:
            raise e

        if "fetch_mRNA_sequences" not in self.execution_times:
            self.execution_times[
                "fetch_mRNA_sequences"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def predict_miRNAs(self, prediction_type: str = "miRDB") -> None:
        logger.info("Started miRNA prediction analysis.")
        self.timer.set_start_time()

        # check the prediction type
        if prediction_type == "miRDB":
            # use the miRDB60predictor to predict miRNAs # TODO make it so that the user submitts the predictior, like metrices
            predictor = miRDB60predictor()
            # iterate through each product and predict miRNAs
            with logging_redirect_tqdm():
                for product in tqdm(self.products, desc="Predict miRNAs"):
                    match_dict = predictor.predict_from_product(
                        product
                    )  # bottleneck operation
                    # if there are matches, add them to the corresponding miRNA objects
                    if match_dict is not None:
                        for miRNA_id, match in match_dict.items():
                            # check if the miRNA already exists in the list of miRNAs
                            for mirna in self.miRNAs:
                                if mirna.id == miRNA_id:
                                    mirna.mRNA_overlaps[product.uniprot_id] = match
                                    break
                            # if the miRNA doesn't exist in the list, create a new miRNA object
                            else:
                                self.miRNAs.append(
                                    miRNA(
                                        miRNA_id,
                                        mRNA_overlaps={product.uniprot_id: match},
                                    )
                                )

        elif prediction_type == "other_type":
            # do something else
            pass
        else:
            # raise an error if the prediction type is invalid
            raise ValueError("Invalid prediction type")

        if "predict_miRNAs" not in self.execution_times:
            self.execution_times["predict_miRNAs"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def change_miRNA_overlap_treshold(
        self, treshold: float, safety: bool = False
    ) -> None:
        """
        Sets the model's 'miRNA_overlap_threshold' to a new 'threshold'. The threshold should be between 0.0 and 1.0.

        WARNING: Changing the miRNA overlap threshold will delete all the calculated previous miRNA scores.

        Parameters:
          - (float) threshold: the new miRNA_overlap_threshold
          - (bool) safety: if False, will ask for the user's confirmation during runtime.
        """
        self.miRNA_overlap_treshold = treshold
        logger.warning(
            "Sorry, but changing the treshold will delete all the calculated miRNA"
            " scores. You will have to calculate them again!"
        )
        if not safety:
            confirmation = input("Are you sure you want to proceed? (y/n)")
            if confirmation.lower() != "y":
                print("Aborting operation.")
                return
        for _miRNA in self.miRNAs:
            _miRNA.scores = {}

    def score_miRNAs(
        self, score_class: List[Metrics], recalculate: bool = False
    ) -> None:
        """
        Performs miRNA scoring on the current ReverseLookup's 'miRNAs' using the input Metrics implementation(s). This function allows the user
        to pass a custom or a pre-defined scoring algorithm, which is of the 'Metrics' type (look in Metrics.py), or a list of scoring algorithms.
        Each miRNA class of the current ReverseLookup instance has a member field 'scores'. For each miRNA instance, score is computed
        and saved to the miRNA's 'scores' dictionary as a mapping between the scoring algorithm's name (eg. "basic_miRNA_score") and the
        corresponding miRNA's float score computed with this scoring algorithm. If multiple scoring algorithms are used, then the miRNA's
        'scores' dictionary will have multiple elements, each a mapping between the scoring algorithm's name and the corresponding score.

        Parameters:
          - score_class: A subclass (implementation) of the Metrics superclass (interface). Current pre-defined Metrics implementations subclasses
                         are 'adv_product_score', 'nterms', 'inhibited_products_id', 'basic_mirna_score'.

                         If 'inhibited_products_id' are used, then the miRNA's 'scoring' field will have a key "inhibited products id", the
                         value at this key will be a list of all the product ids (of all the current GOTerm-associated products, which satisfy
                         the condition that the product's mRNA binding strength > miRNA_overlap_threshold)

                         If 'basic_mirna_score' is used, then [TODO]

          - recalculate: If set to True, will perform score recalculations irrespective of whether a score has already been computed.
                         If set to False, won't perform score recalculations.

        Calling example:
        (1) Construct a ReverseLookup model
        model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

        (2) Create one or more Metrics scoring implementations for the model:
        adv_score = adv_product_score(model)
        nterms_score = nterms(model)

        (3) Call the score_products on the model using the Metrics scoring implementations
        model.score_products([adv_score, nterms_score])
        """
        logger.info("Started miRNA scoring.")
        self.timer.set_start_time()

        if not isinstance(score_class, list):
            score_class = [score_class]

        with logging_redirect_tqdm():
            # iterate over miRNAs using tqdm for progress tracking
            for mirna in tqdm(self.miRNAs, desc="Score miRNAs"):
                # if there is no overlap, skip the miRNA
                if not mirna.mRNA_overlaps:
                    continue
                for _score_class in score_class:
                    if (
                        _score_class.name not in mirna.scores and recalculate is True
                    ):  # if score hasn't been computed, compute it
                        mirna.scores[_score_class.name] = _score_class.metric(mirna)
                    elif _score_class.name not in mirna.scores:
                        mirna.scores[_score_class.name] = _score_class.metric(mirna)

        if "score_miRNAs" not in self.execution_times:
            self.execution_times["score_miRNAs"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    # housekeeping functions

    def get_all_goterms_for_product(self, product: Product | str) -> List[GOTerm]:
        """
        func desc

        Args:
          - (Product) | (str): either a Product object, or a string denoting either a product's UniProtKB id (eg. 'Q8TED9') or a product's
                               gene name (eg. 'AFAP1L1'). A UniProtKB can be input either in the 'UniProtKB:Q8TED9' or the 'Q8TED9' notation.

        Returns:
          - List[GOTerm]: a list of GO Term objects, which are associated with the input Product or product string (UniProtKB id or gene name)
        """
        if isinstance(product, str):
            if ":" in product:
                product = product.split(":")[
                    1
                ]  # if in UniProtKB:xxxx notation, obtain only the last part of the id, eg. 'Q8TED9'
            for prod in self.products:
                if prod.uniprot_id == product:
                    product = prod
                    break
                if prod.genename == product:
                    product = prod
                    break

        goterms_list = []
        for goterm in self.goterms:  # loop over all GO Terms
            if any(
                product_id in goterm.products for product_id in product.id_synonyms
            ):  # a GOTerm has GOTerm.products stored in the full-identifier notation (eg. 'MGI:1201409', 'UniProtKB:Q02763', ...), therefore you need to use product.id_synonyms, which also contains the full-identifier notation
                goterms_list.append(goterm)
        return goterms_list

    def get_all_goterms_for_process(self, process: str) -> List[GOTerm]:
        """
        Loops through all GO Term objects in self.goterms (initialised from input.txt or from load_model at object creation)
        and adds each GO Term instance to a result list, if any of it's processes (goterm.processes) are involved in the parameter 'process'.

        Returns:
          - List[GOTerm]: a list of all GO Term objects, which are associated with the 'process'.

        Example: if process = "diabetes", then it will return a list of all the diabetes-associated GO Terms you specified
        in input.txt file, irrespective of direction (either +, - or 0)
        """
        goterms_list = []
        for goterm in self.goterms:
            if any(proc["process"] == process for proc in goterm.processes):
                goterms_list.append(goterm)
        return goterms_list

    def list_goterms_id(self) -> List[str]:
        """
        Returns a list of all GO term IDs in the GO ontology.
        """
        # Use a list comprehension to extract IDs from the GO terms and return the resulting list
        return [goterm.id for goterm in self.goterms]

    def get_goterm(self, identifier) -> GOTerm:
        """
        Return GOTerm based on any id
        """
        goterm = next(
            obj
            for obj in self.goterms
            if any(
                getattr(obj, attr) == identifier
                for attr in ["id", "name", "description"]
            )
        )
        return goterm

    def get_product(self, identifier) -> Product:
        """
        Return GOTerm based on any id
        """
        product = next(
            obj
            for obj in self.products
            if any(
                getattr(obj, attr) == identifier
                for attr in [
                    "genename",
                    "description",
                    "uniprot_id",
                    "ensg_id",
                    "enst_id",
                    "refseq_nt_id",
                    "mRNA",
                ]
            )
            or any(id == identifier for id in obj.id_synonyms)
        )
        return product

    def save_model(self, filepath: str) -> None:
        data = {}
        data["target_processes"] = self.target_processes
        data["go_categories"] = self.go_categories
        data["filepaths"] = self.datafile_paths
        data["model_settings"] = self.model_settings.to_json()
        data["miRNA_overlap_treshold"] = self.miRNA_overlap_treshold
        data["execution_times"] = self.execution_times
        data["statistically_relevant_products"] = self.statistically_relevant_products

        # save goterms
        for goterm in self.goterms:
            data.setdefault("goterms", []).append(goterm.__dict__)
        # save products
        for product in self.products:
            data.setdefault("products", []).append(product.__dict__)
        # save miRNAs
        for mirna in self.miRNAs:
            data.setdefault("miRNAs", []).append(mirna.__dict__)

        JsonUtil.save_json(data, filepath)

    def compare_to(
        self,
        compare_model: ReverseLookup,
        compare_field: str = "",
        compare_subfields: list = [],
        exclude_http_errors=True,
    ):
        """
        Compares 'compare_field'(s) of this model to the same member fields of 'compare_model'.
        Example: you want to compare if this model has the same GoTerms as the reference 'compare_model': you supply the reference model,
        and set compare_field to "goterms".

        Params:
          - compare_model: a reference ReverseLookup model, against which to compare
          - compare_field: a member field of a ReverseLookup model. Possible options are:
                - 'goterms' - to compare go terms
                - 'products' - to compare products
                - "" (empty) - compare go terms and products in a single function call
                - [TODO]: miRNAs
          - compare_subfields: a list of subfields to compare. For example, if you choose 'goterms' as compare field,
                               you may choose 'name' to compare if the newly server-queried name of a specific go term equals the name of that go term in the reference model.
                - if you choose 'goterms' as compare_field, the options are:
                    - 'name'
                    - 'description'
                    - 'weight'
                    - 'products'
                    note: 'id' (eg. GO:00008286) is not an option, since it is used to carry out comparisons between this model and reference model.
                - if you choose 'products' as compare_field, the options are:
                    - 'id_synonms'
                    - 'description'
                    - 'uniprot_id'
                    - 'ensg_id'
                    - 'enst_id'
                    - 'refseq_nt_id'
                    - 'mRNA'
                    - 'scores_adv-score'
                    - 'scores_nterms'
                    - 'scores_binomial-test'
                    - 'scores_fisher-test'
                    note: 'genename' is not an option, since it is used to carry out comparisons between this model and the reference model.
          - exclude_http_errors: If true, will exclude goterms from comparison, which had known http errors [TODO]

        Returns:

        """

        def compare_json_elements(
            src_json, reference_json, _compare_subfields: list, json_type: str
        ):
            """
            Compares source json ('src_json') to reference json ('reference_json'). All compare_fields are compared.
            'json_type' must be either 'goterms' or 'products'.

            Returns a dictionary of result differences between src_json and reference_json.
            """
            result_diff = {}  # a list of differences
            # if we are looping over go terms, then go terms from src and ref are compared with their 'id' field. If we are doing product comparisons, then products are compared using 'genename'.
            element_identifier = "id" if json_type == "goterms" else "genename"

            count = len(reference_json)
            i = 0
            for ref_element in reference_json:
                logger.debug(f"{i}/{count}")
                i += 1
                # ref_element = json.dumps(ref_element.__dict__) # json conversion, so we can access elements using ['id'] etc.
                current_mismatches = []
                ref_element_id = getattr(ref_element, element_identifier)
                src_element = None
                # find the source element with the same id as reference element
                for src_el in src_json:
                    if getattr(src_el, element_identifier) == ref_element_id:
                        src_element = src_el
                        # src_element = json.dumps(src_el.__dict__)
                        break

                # if no source element is found, note the difference
                if src_element is None:
                    result_diff[ref_element_id] = {
                        "mismatches": ["No source element with same id found."]
                    }
                    continue

                # compare all compare_fields, if any are different between ref_element and src_element, note the difference
                for _compare_subfield in _compare_subfields:
                    # copy ref_element and src_element to preserve original ref_element and src_element for further comparisons. this copy is made, because in case of comparing score fields (eg adv_score), which are nested twice, _ref_element is reassigned the product.scores json "subelement", so inidividual scores, such as adv_score are computed on a one-nested json.
                    if "scores" in _compare_subfield:
                        # _ref_element = ref_element['scores'] #JSON-like approach, this was superseded by the class-based approach
                        # _src_element = src_element['scores']
                        _ref_element = getattr(
                            ref_element, "scores"
                        )  # WARNING: _ref_element is now a JSON
                        _src_element = getattr(
                            src_element, "scores"
                        )  # WARNING: _src_element is now a JSON
                        # convert to class
                        _ref_element_class_placeholder = JsonToClass(str(_ref_element))
                        _src_element_class_placeholder = JsonToClass(str(_src_element))
                        _ref_element = (
                            _ref_element_class_placeholder.object_representation
                        )
                        _src_element = (
                            _src_element_class_placeholder.object_representation
                        )
                        # score-related comparison subfields are sent in the format 'scores_binomial-test'. To convert to the correct one-nested comparison subfield, choose the exact score (the element after _) and replace '-' by '_'
                        # 'scores_adv-score' -> 'adv_score'
                        _compare_subfield = _compare_subfield.split("_")[1].replace(
                            "-", "_"
                        )
                    else:
                        _ref_element = ref_element
                        _src_element = src_element

                    if hasattr(_ref_element, _compare_subfield) and hasattr(
                        _src_element, _compare_subfield
                    ):
                        _ref_element_attr_value = getattr(
                            _ref_element, _compare_subfield
                        )
                        _src_element_attr_value = getattr(
                            _src_element, _compare_subfield
                        )

                        # if ref or src element attr value are classes (namespaces), convert them back to json form; SimpleNamespace is used for type check, since that is the placeholder class used for JSON->class conversion for score jsons
                        # TODO: find out a way how to convert a SimpleNamespace back to JSON. I've tried creating a JsonToClass custom class, which holds the source json, but
                        # _ref_element_attr_value can take up only a specific json segment (eg. when _compare_subfield == fisher_test), _ref_element_attr_value corresponds only to the segment of the json, which is encoded by the "fisher_test".
                        # I cannot obtain such fidelity with access to just source_json.
                        """
                        if isinstance(_ref_element_attr_value, SimpleNamespace):
                            # error: SimpleNamespace is not JSON serializable
                            #_ref_element_attr_value = json.dumps(_ref_element_attr_value.__dict__)
                            #_ref_element_attr_value = json.dumps(vars(_ref_element_attr_value))
                            # test = SimpleNamespaceUtil.simpleNamespace_to_json(_ref_element_attr_value) # TODO: finish this                 
                        if isinstance(_src_element_attr_value, SimpleNamespace):
                            # error: SimpleNamespace is not JSON serializable
                            #_src_element_attr_value = json.dumps(_src_element_attr_value.__dict__)
                            _src_element_attr_value = json.dumps(vars(_src_element_attr_value))
                        """
                        if isinstance(_ref_element_attr_value, list) and isinstance(
                            _src_element_attr_value, list
                        ):
                            """
                            We are dealing with two lists. Check if all elements from _ref_element_attr_value list can be found in _src_element_attr_value
                            """
                            missing_ref_elements_in_src = []
                            missing_src_elements_in_ref = []

                            # check for reference elements in src
                            for ref_e in _ref_element_attr_value:
                                if ref_e not in _src_element_attr_value:
                                    missing_ref_elements_in_src.append(ref_e)

                            # check for src elements in ref
                            for src_e in _src_element_attr_value:
                                if src_e not in _ref_element_attr_value:
                                    missing_src_elements_in_ref.append(src_e)

                            if (
                                missing_ref_elements_in_src != []
                                or missing_src_elements_in_ref != []
                            ):
                                current_mismatches.append(
                                    "Compare field array mismatch for"
                                    f" '{_compare_subfield}'\\n   - missing reference"
                                    " elements in src:"
                                    f" {missing_ref_elements_in_src}\\n    - missing"
                                    " source elements in reference:"
                                    f" {missing_src_elements_in_ref}\\n    - ref ="
                                    f" {_ref_element_attr_value}\\n    - src ="
                                    f" {_src_element_attr_value}"
                                )

                        elif _ref_element_attr_value == _src_element_attr_value:
                            continue  # no mismatch, both are same values

                        else:  # compare field mismatch, values are different
                            current_mismatches.append(
                                f"Compare field mismatch for '{_compare_subfield}': ref"
                                f" = '{_ref_element_attr_value}', src ="
                                f" '{_src_element_attr_value}'"
                            )
                    elif not (
                        hasattr(_ref_element, _compare_subfield)
                        and hasattr(_src_element, _compare_subfield)
                    ):
                        continue  # no mismatch, neither element has this _compare_subfield
                    else:  # one element has _compare_subfield, other doesn't find out which.
                        compare_field_in_ref_element = hasattr(
                            _ref_element, _compare_subfield
                        )
                        compare_field_in_src_element = hasattr(
                            _src_element, _compare_subfield
                        )
                        current_mismatches.append(
                            f"Compare field '{_compare_subfield}' doesn't exist in"
                            " reference or source element. Source element:"
                            f" '{compare_field_in_src_element}', Reference element:"
                            f" '{compare_field_in_ref_element}'"
                        )

                    """ # A JSON-like approach to solving the above class-based approach (which uses hasattr and getattr)
                    if _compare_subfield in _ref_element and _compare_subfield in _src_element: # check if compare_field is equal in ref and src element
                        if _ref_element[_compare_subfield] == _src_element[_compare_subfield]:
                            continue
                        else: # compare field mismatch
                            current_mismatches.append(f"Compare field mismatch for '{_compare_subfield}': ref = {_ref_element[_compare_subfield]} --- src = {_src_element[_compare_subfield]}")
                    elif (_compare_subfield in _ref_element and _compare_subfield not in _src_element) or (_compare_subfield not in _ref_element and _compare_subfield in _src_element): # compare_field is not in ref_element or src_element, find out where
                        compare_field_in_ref_element = _compare_subfield in _ref_element
                        compare_field_in_src_element = _compare_subfield in _src_element
                        current_mismatches.append(f"Compare field '{_compare_subfield}' doesn't exist in reference or source element. Source element: {compare_field_in_src_element}, Reference element: {compare_field_in_ref_element}")
                    """
                    if (
                        current_mismatches != []
                    ):  # append mismatches, if any are found, to result_diff
                        result_diff[ref_element_id] = {"mismatches": current_mismatches}
            # return
            return result_diff

        logger.info("Comparing src json to reference json.")

        allowed_goterms_subfields = ["name", "description", "weight", "products"]
        allowed_products_subfields = [
            "id_synonyms",
            "description",
            "uniprot_id",
            "ensg_id",
            "enst_id",
            "refseq_nt_id",
            "mRNA",
            "scores_adv-score",
            "scores_nterms",
            "scores_binomial-test",
            "scores_fisher-test",
        ]

        if compare_field == "goterms":
            # if all compare_subfields are from allowed_goterms_subfields
            if all(
                compare_subfield
                for compare_subfield in compare_subfields
                if compare_subfield in allowed_goterms_subfields
            ):
                src_json = self.goterms
                ref_json = compare_model.goterms
                _cs = (
                    ["name", "description", "weight", "products"]
                    if compare_subfields == []
                    else compare_subfields
                )
                goterms_diff = compare_json_elements(
                    src_json, ref_json, _compare_subfields=_cs, json_type="goterms"
                )
                return goterms_diff  # the difference in all _compare_subfields across src_json and ref_json goterms
            else:
                logger.error(
                    "Error: one of the supplied compare_subfields"
                    f" ({compare_subfields}) is not allowed for compare field"
                    f" '{compare_field}'. Allowed compare subfields for"
                    f" '{compare_field}' are {allowed_goterms_subfields}"
                )
        elif compare_field == "products":
            # if all compare_subfields are from allowed_products_subfields
            if all(
                compare_subfield
                for compare_subfield in compare_subfields
                if compare_subfield in allowed_products_subfields
            ):
                src_json = self.products
                ref_json = compare_model.products
                # if compare_fields parameter is empty, then use all allowed compare fields, otherwise use parameter
                _cs = (
                    [
                        "id_synonyms",
                        "description",
                        "uniprot_id",
                        "ensg_id",
                        "enst_id",
                        "refseq_nt_id",
                        "mRNA",
                        "scores_adv-score",
                        "scores_nterms",
                        "scores_binomial-test",
                        "scores_fisher-test",
                    ]
                    if compare_subfields == []
                    else compare_subfields
                )
                products_diff = compare_json_elements(
                    src_json, ref_json, _compare_subfields=_cs, json_type="products"
                )
                return products_diff  # the difference in all _compare_subfields across src_json and ref_json products
            else:
                logger.error(
                    "Error: one of the supplied compare_subfields"
                    f" ({compare_subfields}) is not allowed for compare field"
                    f" '{compare_field}'. Allowed compare subfields for"
                    f" '{compare_field}' are {allowed_products_subfields}"
                )
        elif (
            compare_field == ""
        ):  # If compare_field wasn't set, perform comparison on both goterms and products.
            # deduce which compare subfields should be analysed for goterms and which for products
            analysis_goterms_subfields = []  # comparisons will be performed on these
            analysis_products_subfields = []  # comparisons will be performed on these
            for compare_subfield in compare_subfields:
                if compare_subfield in allowed_goterms_subfields:
                    analysis_goterms_subfields.append(compare_subfield)
                elif compare_subfield in allowed_products_subfields:
                    analysis_products_subfields.append(compare_subfield)
                elif (
                    compare_subfield in allowed_goterms_subfields
                    and compare_subfield in allowed_products_subfields
                ):
                    analysis_goterms_subfields.append(compare_subfield)
                    analysis_products_subfields.append(compare_subfield)

            goterms_src_json = self.goterms
            goterms_ref_json = compare_model.goterms
            # use all allowed_goterms_subfields if analysis_goterms_subfields is empty, else use anaylsis_goterms_subfields
            _cs = (
                allowed_goterms_subfields
                if analysis_goterms_subfields == []
                else analysis_goterms_subfields
            )
            goterms_diff = compare_json_elements(
                goterms_src_json,
                goterms_ref_json,
                _compare_subfields=_cs,
                json_type="goterms",
            )

            products_src_json = self.products
            products_ref_json = compare_model.products
            # use all allowed_products_subfields if analysis_products_subfields is empty, else use anaylsis_products_subfields
            _cs = _cs = (
                allowed_products_subfields
                if analysis_products_subfields == []
                else analysis_products_subfields
            )
            products_diff = compare_json_elements(
                products_src_json,
                products_ref_json,
                _compare_subfields=_cs,
                json_type="products",
            )

            # merge both dictionaries
            return {**goterms_diff, **products_diff}

    def perform_statistical_analysis(self, test_name="fisher_test", filepath=""):
        """
        Finds the statistically relevant products, saves them to 'filepath' (if it is provided) and returns a JSON object with the results.

        Parameters:
          - (str) test_name: The name of the statistical test to use for product analysis. It must be either 'fisher_test' (the results of the fisher's test are then used)
          or 'binomial_test' (the results of the binom test are used).
          - (str) filepath: The path to the output file

        Warning: Binomial test scoring is not yet implemented.
        Warning: Products in this model must be scored with the aforementioned statistical tests prior to calling this function.

        Usage example:
            model = ReverseLookup.load_model("diabetes_angio_2/data.json")
            goaf = GOAnnotiationsFile()
            binom_score = binomial_test(model, goaf)
            fisher_score = fisher_exact_test(model, goaf)
            model.score_products([binom_score, fisher_score])
            model.perform_statistical_analysis("fisher")

        Returns a JSON with the following structure (example is also provided to the right):
            {                           {
            PROCESS_PAIR_CODE: [        "diabetes+:angio+": [
                PRODUCT1_DICT               { // product info: id_synonyms, genename, description, ...},
                PRODUCT2_DICT               { // product info: id_synonyms, genename, description, ...},
                ...                         ...
            ],                          ],
            ...                         "diabetes+:obesity+": [...],
                                        "angio+:obesity+": [...]
            }                           }

        The genes (products) for each process pair are sorted according to the sum of the p-values, with products with the lowest pvalues (highest
        statistical probabilities) appearing first in the sorted dictionary.

        TODO: implement binomial score, maybe even adv_score and nterms for backwards compatibility
        """

        def sorting_key(product):
            """
            Sorting key used for the sorting of JSON data based on ascending pvalues.
            Fisher test MUST be calculated for this to work.
            """
            pvalue_sum = 0
            for process in self.target_processes:
                pvalue_process = product["scores"]["fisher_test"][
                    f"{process['process']}{process['direction']}"
                ]["pvalue_corr"]
                pvalue_sum += pvalue_process
            return pvalue_sum

        statistically_relevant_products = (
            []
        )  # a list of lists; each member is [product, "process1_name_direction:process2_name_direction"]
        for product in self.products:
            # given three process: diabetes, angio, obesity, this code iterates through each 2-member combination possible
            #
            # loop iteration \ process      diabetes    angio   obesity
            # it. 0  (i=0,j=1)                  |         |
            # it. 1  (i=0,j=2)                  |                  |
            # it. 2  (i=1,j=2)                            |        |
            # j = 3 -> loop condition not met
            # i = 2 -> loop condition not met (i would start on 'obesity', wouldn't find matching pair with j)
            #
            # Each member pair is used to assess statistically relevant genes, which either positively or
            # negatively regulate both of the processes in the pair.

            for i in range(len(self.target_processes) - 1):
                for j in range(i + 1, len(self.target_processes)):
                    process1 = self.target_processes[i]
                    process2 = self.target_processes[j]
                    pair = [process1, process2]

                    if all(
                        float(
                            product.scores[test_name][
                                f"{process['process']}{process['direction']}"
                            ].get("pvalue_corr", 1)
                        )
                        < self.model_settings.pvalue
                        for process in pair
                    ) and all(
                        float(
                            product.scores[test_name][
                                f"{process['process']}{'+' if process['direction'] == '-' else '-'}"
                            ].get("pvalue_corr", 1)
                        )
                        >= self.model_settings.pvalue
                        for process in pair
                    ):
                        statistically_relevant_products.append(
                            [
                                product,
                                f"{process1['process']}{process1['direction']}:{process2['process']}{process2['direction']}",
                            ]
                        )

        # statistically_relevant_products stores a list of lists, each member list is a Product object bound to a specific pair code (e.g. angio+:diabetes+).
        # statistically_relevant_products_final is a dictionary. It's keys are process pair codes (e.g. angio+:diabetes+), each key holds a list of all statistically relevant products for the process pair
        # (eg. if angio+:diabetes+ it holds all products, which positively regulate both angiogenesis and diabetes)
        process_pairs = (
            []
        )  # each element is a code binding two processes and their direction, eg. angio+:diabetes+
        statistically_relevant_products_final = (
            {}
        )  # dictionary between two processes (eg. angio+:diabetes+) and all statistically relevant products
        for i in range(len(self.target_processes) - 1):
            for j in range(i + 1, len(self.target_processes)):
                process1 = self.target_processes[i]
                process2 = self.target_processes[j]
                pair_code = f"{process1['process']}{process1['direction']}:{process2['process']}{process2['direction']}"
                process_pairs.append(pair_code)
                statistically_relevant_products_final[
                    pair_code
                ] = []  # initialise to empty list

        for element in statistically_relevant_products:
            # each element is a list [product, "process1_name_direction:process2_name_direction"]
            prod = element[0]
            process_pair_code = element[1]
            statistically_relevant_products_final[process_pair_code].append(
                prod.__dict__
            )

        # sort the genes based on the ascending sum of pvalues (lowest pvalues first)
        statistically_relevant_products_final_sorted = {}
        for i in range(len(self.target_processes) - 1):
            for j in range(i + 1, len(self.target_processes)):
                process1 = self.target_processes[i]
                process2 = self.target_processes[j]
                pair_code = f"{process1['process']}{process1['direction']}:{process2['process']}{process2['direction']}"
                statistically_relevant_products_for_process_pair = (
                    statistically_relevant_products_final[pair_code]
                )
                statistically_relevant_products_for_process_pair_sorted = sorted(
                    statistically_relevant_products_for_process_pair,
                    key=lambda gene: sorting_key(gene),
                )
                statistically_relevant_products_final_sorted[
                    pair_code
                ] = statistically_relevant_products_for_process_pair_sorted

        # TODO: save statistical analysis as a part of the model's json and load it up on startup
        self.statistically_relevant_products = (
            statistically_relevant_products_final_sorted
        )
        logger.info(
            "Finished with product statistical analysis. Found"
            f" {len(statistically_relevant_products)} statistically relevant products."
        )
        JsonUtil.save_json(
            data_dictionary=statistically_relevant_products_final_sorted,
            filepath=filepath,
        )
        return statistically_relevant_products_final

    def change_products_member_field(self, member_field_name: str, value):
        """
        This function changes the 'member_field_name' member variable of all Product instances in self.products
        to 'value'.

        Args:
          - (str) member_field_name: The name of the member variable / attribute of a Product instance, the value of which you want to change.
                                     A valid member variable is any member variable of the Product class, such as 'id_synonyms', 'genename', 'had_orthologs_computed' etc
        """
        for product in self.products:
            if hasattr(product, member_field_name):
                setattr(product, member_field_name, value)

    @classmethod
    def load_model(cls, filepath: str) -> "ReverseLookup":
        """
        # TODO: method explanation

        This method also parses GO Term parents for GO Terms, if setting include_indirect_annotations is True.
        """
        data = JsonUtil.load_json(filepath)
        target_processes = data["target_processes"]
        miRNA_overlap_treshold = data["miRNA_overlap_treshold"]

        execution_times = {}
        if "execution_times" in data:
            execution_times = data["execution_times"]

        if "statistically_relevant_products" in data:
            statistically_relevant_products = data["statistically_relevant_products"]
        else:
            statistically_relevant_products = {}

        if "go_categories" in data:
            go_categories = data["go_categories"]
        else:
            go_categories = [
                "biological_process",
                "molecular_activity",
                "cellular_component",
            ]

        if "model_settings" in data:
            settings = ModelSettings.from_json(data["model_settings"])
        else:
            settings = ModelSettings()

        if "filepaths" in data:
            datafile_paths = data["filepaths"]
        else:
            datafile_paths = {}

        goterms = []
        for goterm_dict in data["goterms"]:
            goterms.append(GOTerm.from_dict(goterm_dict))

        obo_parser = None
        if settings.include_indirect_annotations is True:
            if datafile_paths != {} and "go_obo_filepath" in datafile_paths:
                if (
                    datafile_paths["go_obo_filepath"] is not None
                    and datafile_paths["go_obo_filepath"] != ""
                ):
                    obo_parser = OboParser(datafile_paths["go_obo_filepath"])
            else:
                obo_parser = OboParser()
            for goterm in goterms:
                assert isinstance(goterm, GOTerm)
                if goterm.parent_term_ids == [] or goterm.parent_term_ids is None:
                    goterm_obo = obo_parser.all_goterms[
                        goterm.id
                    ]  # obo representation of this goterm
                    goterm.update(
                        goterm_obo
                    )  # update current goterm with information from .obo file

                    goterm_parent_ids = obo_parser.get_parent_terms(
                        goterm.id
                    )  # calculate parent term ids for this goterm
                    goterm.parent_term_ids = goterm_parent_ids  # update parent term ids

        products = []
        for product_dict in data.get("products", []):
            products.append(Product.from_dict(product_dict))

        miRNAs = []
        for miRNAs_dict in data.get("miRNAs", []):
            miRNAs.append(miRNA.from_dict(miRNAs_dict))

        return cls(
            goterms,
            target_processes,
            products,
            miRNAs,
            miRNA_overlap_treshold,
            execution_times=execution_times,
            statistically_relevant_products=statistically_relevant_products,
            go_categories=go_categories,
            model_settings=settings,
            obo_parser=obo_parser,
            datafile_paths=datafile_paths,
        )

    @classmethod
    def from_input_file(cls, filepath: str) -> "ReverseLookup":
        """
        Creates a ReverseLookup object from a text file.

        Args:
            filepath (str): The path to the input text file.

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        # Define constants used in parsing the file
        LINE_ELEMENT_DELIMITER = "\t"  # Data is tab separated
        COMMENT_DELIMITER = "#"  # Character used to denote a comment
        LOGIC_LINE_DELIMITER = "###" # Special set of characters to denote a "logic line"

        target_processes = []
        go_categories = []
        go_terms = []
        settings = ModelSettings()
        datafile_paths = {} # a dictionary between data file types and their filepaths {'go_obo_filepath': "GO_OBO_FILEPATH", {...}, ...}

        def process_comment(line):
            """
            Processes a comment in the line: returns the part of the line before the comment. The input file should be structured to contain
            three sections - 'settings', 'processes' and 'GO_terms', annotated using the LOGIC_LINE_DELIMITER.

            For the construction of input.txt, please refer to the Readme file. [TODO]

            Parameters:
            - line: the line whose comment to process
            """
            if LOGIC_LINE_DELIMITER in line:
                # Logic lines should be marked with "###" at the start. For a logic line, the returned result is line without the line_keep_delimiter
                return line.replace(LOGIC_LINE_DELIMITER, "")

            if COMMENT_DELIMITER in line:
                return line.split(COMMENT_DELIMITER)[0]
            else:
                return line

        filepath_readlines = 0
        with open(filepath, 'r') as f:
            filepath_readlines = len(f.readlines())

        def process_file(filepath: str):
            with open(filepath, "r") as read_content:
                read_lines = read_content.read().splitlines()[2:]  # skip first 2 lines
                section = ""  # what is the current section i am reading
                for line in read_lines:
                    line = process_comment(line)
                    line = line.strip()
                    if line == "":
                        continue
                    if "settings" in line:
                        section = "settings"
                        continue
                    elif "filepaths" in line:
                        section = "filepaths"
                        continue
                    elif "processes" in line:
                        section = "process"
                        continue
                    elif "categories" in line:
                        section = "categories"
                        continue
                    elif "GO_terms" in line:
                        section = "GO"
                        continue
                    if section == "settings":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        setting_name = chunks[0]
                        setting_value = chunks[1]  # is string now
                        if setting_value == "True" or setting_value == "true":
                            setting_value = True
                        if setting_value == "False" or setting_value == "false":
                            setting_value = False
                        if setting_name == "pvalue":
                            setting_value = float(chunks[1])
                        settings.set_setting(
                            setting_name=setting_name, setting_value=setting_value
                        )
                    elif section == "filepaths":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        datafile_paths[chunks[0]] = chunks[1]
                    elif section == "process":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        target_processes.append(
                            {"process": chunks[0], "direction": chunks[1]}
                        )
                    elif section == "categories":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        category = chunks[0]
                        category_value = chunks[1]
                        if category_value == "True":
                            go_categories.append(category)
                    elif section == "GO":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        if len(chunks) == 5:
                            d = {
                                "id": chunks[0],
                                "processes": {
                                    "process": chunks[1],
                                    "direction": chunks[2],
                                },
                                "weight": int(chunks[3]),
                                "name": chunks[4],
                            }
                        else:
                            d = {
                                "id": chunks[0],
                                "processes": {
                                    "process": chunks[1],
                                    "direction": chunks[2],
                                },
                                "weight": int(chunks[3]),
                            }
                        if not any(
                            d["id"] == goterm.id for goterm in go_terms
                        ):  # TODO: check this !!!!!
                            go_terms.append(GOTerm.from_dict(d))
                        else:  # TODO: check this !!!!!
                            next(
                                goterm for goterm in go_terms if d["id"] == goterm.id
                            ).add_process(
                                {"process": chunks[1], "direction": chunks[2]}
                            )

        try:
            process_file(filepath)
        except OSError:
            logger.error(f"ERROR while opening filepath {filepath}")
            return

        # PROCESS INPUT FILE
        """
        if not os.path.isabs(filepath): # this process with traceback.extract_stack works correctly on mac, but not on windows.
            current_dir = os.path.dirname(os.path.abspath(traceback.extract_stack()[0].filename))
            mac_filepath = os.path.join(current_dir, filepath) 
        try:
            os.makedirs(os.path.dirname(mac_filepath), exist_ok=True) # this approach works on a mac computer
            process_file(mac_filepath)
        except OSError:
            # # first pass is allowed, on Windows 10 this tries to create a file at 
            # 'C:\\Program Files\\Python310\\lib\\diabetes_angio_1/general.txt'
            # which raises a permission error.
            # fallback if the above fails
            try:
                win_filepath = FileUtil.find_win_abs_filepath(filepath)
                os.makedirs(os.path.dirname(win_filepath), exist_ok=True)
                process_file(win_filepath)
            except OSError:
                logger.error(f"ERROR while opening win filepath {win_filepath}")
                return
        """
        obo_parser = None
        if settings.include_indirect_annotations:
            if datafile_paths != {} and "go_obo_filepath" in datafile_paths:
                if (
                    datafile_paths["go_obo_filepath"] is not None
                    and datafile_paths["go_obo_filepath"] != ""
                ):
                    obo_parser = OboParser(datafile_paths["go_obo_filepath"])
            else:
                obo_parser = OboParser()
            logger.info(
                "Starting OboParser to find all GO Term parents and children using"
                f" data file {obo_parser.filepath}"
            )

            # update goterms to include all parents and children
            with logging_redirect_tqdm():
                for goterm in tqdm(go_terms, desc="Compute indirect nodes"):
                    assert isinstance(goterm, GOTerm)
                    if goterm.parent_term_ids == [] or goterm.parent_term_ids is None:
                        goterm_obo = obo_parser.all_goterms[
                            goterm.id
                        ]  # obo representation of this goterm
                        goterm.update(
                            goterm_obo
                        )  # update current goterm with information from .obo file

                        goterm_parent_ids = obo_parser.get_parent_terms(
                            goterm.id
                        )  # calculate parent term ids for this goterm
                        goterm_children_ids = obo_parser.get_child_terms(
                            goterm.id
                        )  # calculdate child term ids for this goterm
                        goterm.parent_term_ids = (
                            goterm_parent_ids  # update parent term ids
                        )
                        goterm.child_term_ids = (
                            goterm_children_ids  # update child term ids
                        )
            logger.info("Indirect annotations have been computed.")

            """      
            for goterm in go_terms:
                assert isinstance(goterm, GOTerm)
                if goterm.parent_term_ids == [] or goterm.parent_term_ids == None:
                    goterm_obo = obo_parser.all_goterms[goterm.id] # obo representation of this goterm
                    goterm.update(goterm_obo) # update current goterm with information from .obo file

                    goterm_parent_ids = obo_parser.get_parent_terms(goterm.id) # calculate parent term ids for this goterm
                    goterm_children_ids = obo_parser.get_child_terms(goterm.id) # calculdate child term ids for this goterm
                    goterm.parent_term_ids = goterm_parent_ids # update parent term ids
                    goterm.child_term_ids = goterm_children_ids # update child term ids
            """

        logger.info("Creating model from input file with:")
        logger.info(f"  - input file filepath: {filepath}")
        logger.info(f"  - input file line count: {filepath_readlines}")
        logger.info(f"  - count GO Terms: {len(go_terms)} ")
        logger.info(f"  - target_processes: {target_processes}")
        logger.info(f"  - GO categories: {go_categories}")
        logger.info(f"  - model settings: {settings.to_json()}")
        logger.info(f"  - obo_parser: {obo_parser}")
        logger.info(f"  - datafile_paths: {datafile_paths}")
        return cls(
            go_terms,
            target_processes=target_processes,
            go_categories=go_categories,
            model_settings=settings,
            obo_parser=obo_parser,
            datafile_paths=datafile_paths,
        )

    @classmethod
    def from_dict(cls, data: Dict[str, List[Dict]]) -> "ReverseLookup":
        """
        Creates a ReverseLookup object from a dictionary.

        Args:
            data (dict): A dictionary containing lists of GOTerm and target_processes.

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        goterms = [GOTerm.from_dict(d) for d in data["goterms"]]
        target_processes = data["target_processes"]
        if "go_categories" in data:
            go_categories = data["go_categories"]
        else:
            go_categories = [
                "biological_process",
                "molecular_activity",
                "cellular_component",
            ]
        if "model_settings" in data:
            settings = ModelSettings.from_json(data["model_settings"])
        else:
            settings = ModelSettings()
        if "filepaths" in data:
            datafile_paths = data["filepaths"]
        else:
            datafile_paths = []

        logger.info("Model creation from input file complete.")
        return cls(
            goterms,
            target_processes,
            go_categories=go_categories,
            model_settings=settings,
            datafile_paths=datafile_paths,
        )

    def _debug_shorten_GO_terms(self, count):
        """
        Shortens the amount of GO terms to the specified 'count', for debugging purposes.
        """
        if count < len(self.goterms):
            self.goterms = self.goterms[0:count]

from __future__ import annotations
import asyncio
import aiohttp
from contextlib import asynccontextmanager
from typing import List, Dict
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import json
import os
from statsmodels.stats.multitest import multipletests

from .core.ModelSettings import ModelSettings, OrganismInfo
from .core.ModelStats import ModelStats
from .core.GOTerm import GOTerm
from .core.Product import Product
from .core.Metrics import Metrics, basic_mirna_score, miRDB60predictor
from .core.miRNA import miRNA
from .web_apis.GOApi import GOApi
from .web_apis.EnsemblApi import EnsemblApi
from .web_apis.UniprotApi import UniProtApi
from .web_apis.gProfilerApi import gProfiler
from .parse.GOAFParser import GOAnnotationsFile
from .parse.OrthologParsers import HumanOrthologFinder
from .parse.OBOParser import OboParser
from .util.JsonUtil import JsonUtil, JsonToClass
from .util.Timer import Timer
from .util.WebsiteParser import WebsiteParser
from .util.DictUtil import DictUtil
from .util.FileUtil import FileUtil
from .util.CacheUtil import Cacher
from .util.ApiUtil import EnsemblUtil

import logging

# from logging import config
# config.fileConfig("logging_config.py")
logger = logging.getLogger(__name__)


class TargetSOI:
    def __init__(self, name: str, direction: str) -> None:
        """
        A class representing a target SOI. NOT USED CURRENTLY

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
        target_SOIs: List[Dict[str, str]],
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
        input_filepath:str = None,
        GO_api_version:str = None,
        OBO_version_info:dict = None,
        defined_SOIs: List[Dict[str,str]] = None
    ):
        """
        A class representing a reverse lookup for gene products and their associated Gene Ontology terms.

        Args:
            goterms (set): A set of GOTerm objects.
            target_SOIs (list): A list of dictionaries containing SOI names and directions.
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
            defined_SOIs (list): A list of all target SOIs and, if defined, their reverse SOIs
        """
        self.goterms = goterms
        self.products = products
        self.target_SOIs = target_SOIs
        self.defined_SOIs = defined_SOIs
        self.miRNAs = miRNAs
        self.miRNA_overlap_treshold = miRNA_overlap_treshold
        self.model_settings = model_settings
        self.execution_times = execution_times  # dict of execution times, logs of runtime for functions
        self.timer = Timer()
        self.GO_api_version = GO_api_version
        
        if input_filepath is not None:
            self.input_filepath = input_filepath
            if self.model_settings.destination_dir is None:
                self.model_settings.destination_dir = os.path.dirname(os.path.realpath(input_filepath))
        
        """
        if destination_dir is not None:
            self.destination_dir = destination_dir
        
        if destination_dir is None:
            appdata_dir = os.path.join(os.getenv("APPDATA"), "goreverselookup")
            appdata_dir.replace("\\", "/")
            self.destination_dir = appdata_dir
        """
                
        # placeholder to populate after perform_statistical_analysis is called
        self.statistically_relevant_products = statistically_relevant_products

        # GO categories - determines which categories of GO terms are chosen during the Fisher score computation (either from the GO Annotations File or from a gene id to GO Terms API query)
        self.go_categories = go_categories

        self.goaf = GOAnnotationsFile(
            filepath=self.model_settings.get_datafile_path("goa_human"), 
            go_categories=go_categories, 
            valid_evidence_codes=model_settings.valid_evidence_codes, 
            evidence_codes_to_ecoids=model_settings.evidence_codes_to_ecoids
        )
        
        if self.goaf is None:
            logger.warning("MODEL COULD NOT CREATE A GO ANNOTATIONS FILE!")
            logger.warning(f"  - goa_human = {self.model_settings.get_datafile_path('goa_human')}")

        self.go_api = GOApi() # this enables us to use go_api inside Metrics.py, as importing GOApi inside Metrics.py creates circular imports.

        if self.GO_api_version is None:
            self.GO_api_version = self.go_api.get_GO_version()

        if obo_parser is not None:
            self.obo_parser = obo_parser
        else:
            # obo parser is none
            if (
                self.model_settings.datafile_paths.get("go_obo") is not None
                and self.model_settings.datafile_paths["go_obo"]["local_filepath"] != ""
            ):
                self.obo_parser = OboParser(obo_filepath=self.model_settings.get_datafile_path('go_obo'), obo_download_url=self.model_settings.get_datafile_url('go_obo'))
            else:
                self.obo_parser = OboParser()
        
        if OBO_version_info is None:
            OBO_version_info = {
                'format-version': self.obo_parser.format_version,
                'data-version': self.obo_parser.data_version,
                'ontology': self.obo_parser.ontology
            }
        self.OBO_version_info = OBO_version_info
        
        ModelStats.goterm_count = len(self.goterms)
        ModelStats.product_count = len(self.products)

        self.total_goterms = len(self.goterms)
        self.total_products = len(self.products)

        WebsiteParser.init()


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

        if "fetch_all_go_term_names_descriptions" not in self.execution_times:  # to prevent overwriting on additional runs of the same model name
            self.execution_times["fetch_all_go_term_names_descriptions"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_all_go_term_names_descriptions_async(
        self, 
        api: GOApi, 
        req_delay=0.1, 
        max_connections=50
    ):
        """
        Call fetch_all_go_term_names_descriptions with run_async == True to run this code.
        """
        connector = aiohttp.TCPConnector(limit=max_connections, limit_per_host=max_connections) 
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
        request_params={"rows": 10000000},
        max_connections=60
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
            source = self.goaf

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
                        model_settings=self.model_settings,
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
                        goterm.fetch_products(source, model_settings=self.model_settings)

        if "fetch_all_go_term_products" not in self.execution_times:
            self.execution_times["fetch_all_go_term_products"] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    async def _fetch_all_go_term_products_async_v1(
        self, 
        recalculate: bool = False, 
        delay: float = 0.0
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
        self, max_connections=100, request_params={"rows": 1000000}, req_delay=0.5
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

        connector = aiohttp.TCPConnector(limit=max_connections, limit_per_host=max_connections)  # default is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            for goterm in self.goterms:
                url = api.get_products(
                    goterm.id, get_url_only=True, request_params=request_params
                )
                await asyncio.sleep(req_delay)  # request delay
                response = await session.get(url)

                if response.status != 200:  # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                    logger.warning(f"HTTP Error when parsing {goterm.id}. Response status = {response.status}")
                    goterm.http_error_codes["products"] = (f"HTTP Error: status = {response.status}, reason = {response.reason}")

                # data = await response.json()
                response_content = await response.read()
                data = json.loads(response_content)

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

    async def _fetch_all_goterm_products_async_v3(
        self,
        model_settings:ModelSettings,
        max_connections=100,
        request_params={"rows": 10000000},
        req_delay=0.5,
        recalculate: bool = False,
    ):
        """
        In comparison to (GOApi)._fetch_all_go_term_products_async, this function doesn't overload the server and cause the server to block our requests.
        In comparison to the v2 version of this function (inside GOApi), v3 uses asyncio.gather, which speeds up the async requests.
        """
        connector = aiohttp.TCPConnector(limit=max_connections, limit_per_host=max_connections)  # default is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = []
            for goterm in self.goterms:
                if goterm.products == [] or recalculate is True:
                    task = goterm.fetch_products_async_v3(
                        session, model_settings=model_settings, request_params=request_params, req_delay=req_delay
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
        logger.info(f"Creating products from GO Terms. Num goterms = {len(self.goterms)}")
        self.timer.set_start_time()

        # Create an empty set to store unique products
        products_set = set()

        # determine target taxon
        target_taxon = None
        if self.model_settings is not None:
            target_taxon = self.model_settings.target_organism.ncbi_id_full

        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            if term is None:
                continue
            if term.products is None:
                continue
            for product in term.products:
                if product not in products_set:
                    products_set.update(product)
                    if ":" in product:
                        product_object = Product.from_dict({"id_synonyms": [product], "taxon": term.products_taxa_dict[product], "target_taxon": target_taxon})
                    else:
                        product_object = Product.from_dict({"id_synonyms": [product], "genename": product, "taxon": term.products_taxa_dict[product], "target_taxon": target_taxon})
                    self.products.append(product_object)

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

        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            products_set.update(term.products)

        # Iterate over each product in the products_set and create a new Product object from the product ID using the
        # Product.from_dict() classmethod. Add the resulting Product objects to the ReverseLookup object's products list.
        i = 0
        for product in products_set:  # here, each product is a product id, eg. 'MGI:1343124'
            if check_exists(product) is True:
                continue
            if ":" in product:
                self.products.append(Product.from_dict({"id_synonyms": [product]}))
            else:
                self.products.append(Product.from_dict({"id_synonyms": [product], "genename": product}))
            i += 1

        ModelStats.product_count = len(self.products)
        logger.info(f"Created {i} Product objects from GOTerm object definitions")
        logger.info(f"Total Product objects in ReverseLookup model: {len(self.products)}")

        if "create_products_from_goterms" not in self.execution_times:
            self.execution_times[
                "create_products_from_goterms"
            ] = self.timer.get_elapsed_formatted()
        self.timer.print_elapsed_time()

    def fetch_uniprot_product_ensg_ids(self):
        """
        Fetches the ENSG identifiers for the gene products, which have a UniProtKB identifier.
        This operation relies on a bulk request to the UniProtKB servers to perform the id mapping between UniProtKB
        and Ensembl. It is advisable to perform it before attempting the ortholog search.
        """
        uniprot_api = UniProtApi()

        if self.products == [] or self.products is None:
            if len(self.goterms) > 0:
                self.create_products_from_goterms()
            if self.products == [] or self.products is None:
                raise Exception("Model has no defined gene products!")
        
        uniprot_product_ids = []
        for product in self.products:
            assert isinstance(product, Product)
            if product.uniprot_id == None:
                continue
            if "UniProtKB" in product.uniprot_id:
                uniprot_product_ids.append(product.uniprot_id)
        
        logger.info(f"Performing batch UniProtKB id to Ensembl gene id mapping for {len(uniprot_product_ids)} genes.")
        uniprot_to_ensembl_idmap = uniprot_api.idmapping_ensembl_batch(uniprot_ids=uniprot_product_ids)
        if uniprot_to_ensembl_idmap is not None:
            successful_conversions = uniprot_to_ensembl_idmap['results']
            failed_conversions = uniprot_to_ensembl_idmap['failedIds']
            logger.info(f"UniProtKB->Ensembl id mapping results: successful = {len(successful_conversions)}, failed = {len(failed_conversions)}")
        else:
            logger.error(f"UniProtKB->Ensembl id mapping wasn't performed due to an error.")
            return

        target_organism_stable_id_prefix = EnsemblUtil.taxon_to_ensembl_stable_id_prefix(self.model_settings.target_organism.ncbi_id_full)

        # update self.produts
        for idmap in successful_conversions:
            uniprot_id = idmap['from']
            ensembl_id = idmap['to']
            if "." in ensembl_id: # ENSGxxxx.1 -> ENSGxxxx
                ensembl_id = ensembl_id.split(".")[0]
            ens_stable_id_prefix = EnsemblUtil.split_ensembl_id(ensembl_id).get('stable_id_prefix')
            for product in self.products:
                if product.uniprot_id == None:
                    continue
                if uniprot_id in product.uniprot_id:
                    if ens_stable_id_prefix == target_organism_stable_id_prefix:
                        product.ensg_id = ensembl_id
                    else:
                        product.id_synonyms.append(ensembl_id)
    
    def products_perform_idmapping(
        self,
        from_databases:list = ["UniProtKB", "ZFIN", "MGI", "RGD", "Xenbase"],
        to_db = "UniProtKB"
    ):
        """
        First, attempts to map all non-UniProtKB gene product ids (from self.products) to respective UniProtKB gene ids (by default, since 'to_db' is set to UniProtKB).
        In the second run, attempts to map all UniProtKB gene ids to Ensembl gene ids via a call to (ReverseLookup).fetch_uniprot_product_ensg_ids()

        The allowed from_databases and to_db values can be found using the link: https://rest.uniprot.org/configure/idmapping/fields?content-type=application/json
        """
        uniprot_api = UniProtApi()

        if self.products == [] or self.products is None:
            if len(self.goterms) > 0:
                self.create_products_from_goterms()
            if self.products == [] or self.products is None:
                raise Exception("Model has no defined gene products!")

        ### Run 1: map all non-uniprotkb ids to uniprotkb ids
        non_uniprot_dbs = [db for db in from_databases if "UniProtKB" not in db] # select all dbs besides uniprot
        non_uniprot_dbs_product_dict = {} # mapping {'ZFIN': [list_of_zfin_product_ids], ...}
        non_uniprot_dbs_product_dict_new = {}
        
        non_uniprot_product_dict_cached_successful = list()
        non_uniprot_product_dict_cached_failed = list()
        for db in non_uniprot_dbs: # initialise with empty values
            non_uniprot_dbs_product_dict[db] = set()
            non_uniprot_dbs_product_dict_new[db] = list()
        
        for product in self.products: # populate non_uniprot_dbs_product_dict
            if "UniProtKB" in product.id_synonyms[0]:
                continue # don't perform uniprotkb-uniprotkb mappings
            p_id = product.id_synonyms[0] # ZFIN:XXXX
            p_id_db = p_id.split(":")[0] # ZFIN
            try:
                non_uniprot_dbs_product_dict[p_id_db].update([p_id])
            except KeyError:
                logger.debug(f"KeyError when trying to search for {p_id_db} in non_uniprot_dbs_product_dict. This is benign.")
        
        # convert sets to lists
        for db, db_set in non_uniprot_dbs_product_dict.items():
            non_uniprot_dbs_product_dict[db] = list(db_set)
        
        # exclude cached mappings
        for db, db_list in non_uniprot_dbs_product_dict.items():
            start_count = len(db_list)
            new_ids = []
            
            for id in db_list:
                id = id.split(":")[1]
                data_key_successful = f"product_idmapping[id={id},from_db={db},to_db={to_db}]"
                data_key_failed = f"product_idmapping[id={id},failed-mapping]"
                prev_data_successful = Cacher.get_data(data_location="uniprot", data_key=data_key_successful)
                prev_data_failed = Cacher.get_data(data_location="uniprot", data_key=data_key_failed)
                
                is_new = True
                if prev_data_successful is not None:
                    non_uniprot_product_dict_cached_successful.append(prev_data_successful)
                    is_new = False
                if prev_data_failed is not None:
                    non_uniprot_product_dict_cached_failed.append(prev_data_failed.get("id"))
                    is_new = False
                
                if is_new:
                    new_ids.append(id)
                    
            non_uniprot_dbs_product_dict[db] = new_ids
            end_count = len(non_uniprot_dbs_product_dict[db])
            logger.info(f"Cached {start_count-end_count} ids for idmapping for database {db}. Start count: {start_count}, end_count (these will be queried): {end_count}")
        
        # perform batch id mapping
        for db, db_ids_list in non_uniprot_dbs_product_dict.items():
            if db_ids_list == []:
                logger.info(f"Skipping {db}, as there are no associated genes.")
                continue
            idmappings = uniprot_api.idmapping_batch(
                ids = db_ids_list,
                from_db=db,
                to_db=to_db
            )
            if idmappings is None:
                continue
            successful_mappings = idmappings['results']
            failed_mappings = idmappings['failedIds']
            logger.info(f"{db} id mapping conversion to UniProtKB: {len(successful_mappings)} successful mappings, {len(failed_mappings)} failed mappings.")
            for idmap in successful_mappings:
                initial_id = idmap['from']
                uniprot_id = idmap['to']['primaryAccession']
                uniprot_id_full = f"UniProtKB:{uniprot_id}"
                
                data_key = f"product_idmapping[id={initial_id},from_db={db},to_db={to_db}]"
                Cacher.store_data(
                    data_location="uniprot",
                    data_key=data_key,
                    data_value={
                        'from_id': initial_id,
                        'to_id': uniprot_id,
                        'from_db': db,
                        'to_db': to_db,
                        'status': "successful"
                    }
                )
                
                for product in self.products:
                    if len(product.id_synonyms) == 0:
                        continue
                    if initial_id in product.id_synonyms[0] and db in product.id_synonyms[0]:
                        product.uniprot_id = uniprot_id_full
            
            # store failed ids in cache
            for id in failed_mappings:
                data_key = f"product_idmapping[id={id},failed-mapping]"
                Cacher.store_data(
                    data_location="uniprot", 
                    data_key=data_key, 
                    data_value= {
                        'id': id,
                        'status': "failed"
                    }
                )
        
        # also process cached ids
        for element in non_uniprot_product_dict_cached_successful:
            initial_id = element['from_id']
            mapped_id = element['to_id']
            to_db = element['to_db']
            mapped_id_full = f"{to_db}:{mapped_id}"
            for product in self.products:
                if len(product.id_synonyms) == 0:
                        continue
                for syn in product.id_synonyms: # check all synonyms for a match
                    if initial_id in syn and db in syn:
                        if "UniProtKB" in mapped_id_full:
                            product.uniprot_id = mapped_id_full         
        
        ### Run 2: attempt Ensembl gene mappings
        self.fetch_uniprot_product_ensg_ids()

    def fetch_orthologs_products_batch_gOrth(self, target_taxon_number="9606"):
        """
        Fetches a whole batch of gene product orthologs using gOrth (https://biit.cs.ut.ee/gprofiler/orth).

        If target organism was defined in the input.txt file to the program, then 'target_taxon_number' will be overriden
        by the taxon number specified for the target organism.

        Note: This function automatically uses the ENTREZGENE-ACC gOrth namespace (list of all namespaces: https://biit.cs.ut.ee/gprofiler/page/namespaces-list).
        In general, the ENTREZGENE_ACC gOrth namespace works correctly with the following tested identifier types:
          - UniProtKB
          - MGI

        WARNING: The following identifiers are known NOT to work with the ENTREZGENE-ACC namespace:
          - RGD
        For these identifiers, it is highly advisable to call (ReverseLookup).products_perform_idmapping() prior to calling
        fetch_orthologs_products_batch_gOrth, so that RGD ids will be mapped to their respective UniProt and Ensembl ids, which work
        with the ENTREZGENE-ACC namespace.
        """
        if self.model_settings.ortholog_organisms is None or self.model_settings.ortholog_organisms == []:
            logger.info(f"No orthologous organisms specified. Skipping batch gOrth query.")
            return
        
        # determine target_taxon_number if target organism is defined in ModelSettings
        if self.model_settings is not None:
            target_taxon_number = f"{self.model_settings.target_organism.ncbi_id}"

        # divide products into distinct groups
        taxa_ids = set()
        taxa_ids.add(f"{self.model_settings.target_organism.ncbi_id}")
        if self.model_settings.ortholog_organisms is not None:
            for label,ortholog_organism in self.model_settings.ortholog_organisms.items():
                taxa_ids.add(f"{ortholog_organism.ncbi_id}")
        taxa_ids = list(taxa_ids)
        
        # initialise empty dict; format: {'taxon_id_number': [LIST OF ASSOCIATED TAXON IDs]}
        ortholog_query_dict = {}
        for taxon_number in taxa_ids:
            ortholog_query_dict[taxon_number] = []
        
        # iterate over products, construct ortholog_query_dict with id synonyms, uniprot ids and ensembl ids
        for product in self.products:  
            if product.taxon is None:
                continue
            
            product_taxon_number = product.taxon.split(":")[1] if ":" in product.taxon else product.taxon

            if int(product_taxon_number) == self.model_settings.target_organism.ncbi_id: # don't query orthologs if product taxon is the same as target organism taxon
                continue

            product_ids_short = [] # mustn't be UniProtKB:xxxx but just 'xxxx'
            if len(product.id_synonyms) > 0: # search id synonyms
                id_syn = product.id_synonyms[0]
                if ":" in id_syn:
                    p_id_short = id_syn.split(":")[1]
                    product_ids_short.append(p_id_short)
            if product.uniprot_id is not None: # add uniprot query if uniprot for id synonym has been precomputed
                product_ids_short.append(product.uniprot_id.split(":")[1]) if ":" in product.uniprot_id else product_ids_short.append(f"{product.uniprot_id}")     
            if product.ensg_id is not None: # add ensembl query if ensembl id for id synonym has been precomputed
                # ENSRNOG0000145.2 -> ENSRNOG0000145
                product_ids_short.append(product.ensg_id.split(".")[0]) if "." in product.ensg_id else product_ids_short.append(f"{product.ensg_id}")
            if product_ids_short == []: 
                # continue if no id has been found
                continue
            
            for p_id in product_ids_short:
                ortholog_query_dict[product_taxon_number].append(p_id)

        num_queried_orthologs_definitive = 0
        num_queried_orthologs_indefinitive = 0
        num_no_orthologs = 0

        total_ortholog_query_count = 0
        max_single_ortholog_query_count = 0
        for taxon_id, input_ids_to_query in ortholog_query_dict.items():
            total_ortholog_query_count += len(input_ids_to_query)
            max_single_ortholog_query_count = max(max_single_ortholog_query_count, len(input_ids_to_query))

        logger.info(f"gOrth query begin: querying a total of {total_ortholog_query_count} ids, maximum query count per request: {max_single_ortholog_query_count}.")
        # perform search and modify self.products with the search results
        gprofiler = gProfiler()
        logger.info(f"Processing gOrth ortholog query dicts")
        timer = Timer()
        for taxon, ids in ortholog_query_dict.items():
            logger.info(f"  - taxon {taxon}: {len(ids)} items")
            timer.set_start_time()

            ortholog_query_results = gprofiler.find_orthologs(source_ids=ids, source_taxon=taxon, target_taxon=target_taxon_number) #all ortholog query results for this taxon
            if ortholog_query_results is None:
                # None can be returned from gprofiler.find_orthologs if source_ids is []
                continue

            # process results for this taxon - this takes TOO long!
            for id, ortholog_results in ortholog_query_results.items():
                p = self.get_product(id, identifier_type="id_synonyms")
                if p is None:
                    continue

                if ortholog_results == [] and p.gorth_ortholog_exists != True: # p.gorth_ortholog_exists != True is set if this product had already been determined to have an existing ortholog. For example, ENSRNOGxxxx determines a valid gOrth ortholog, sets it, but then RGD:xxxx (pointing to the same gene) finds no gOrth ortholog -> RGD results in this case shouldn't disturb the successful ensembl gOrth ortholog query.
                    p.gorth_ortholog_exists = False
                    p.gorth_ortholog_status = "none"

                if len(ortholog_results) > 1 and p.gorth_ortholog_status != "definitive": # if a valid gOrth ortholog has been found, overwrite any prior queries that didn't find an ortholog for this product. Don't overwrite and prior "definitive"-ly queried gOrth orthologs!
                    # more than 1 ortholog ids found = "indefinitive ortholog"
                    p.gorth_ortholog_exists = True
                    p.gorth_ortholog_status = "indefinitive"

                    # determine autoselect first indefinitive ortholog (based on model settings)
                    autoselect_first_among_multiple_ensg_orthologs = False
                    if self.model_settings is not None:
                        if self.model_settings.gorth_ortholog_fetch_for_indefinitive_orthologs == False:
                            # if ortholog fetch for indefinitive orthologs is False -> we don't want to query orthologs in the
                            # regular ortho query pipeline -> set autoselect to true -> select 1st returned ENSG by gOrth -> no regular ortho pipeline, as
                            # ENSG is already found!
                            autoselect_first_among_multiple_ensg_orthologs = True

                    if autoselect_first_among_multiple_ensg_orthologs:
                        if p.ensg_id is not None: # move existing ensembl id among id synonyms
                            p.id_synonyms.append(p.ensg_id)
                        p.ensg_id = ortholog_results[0]
                        # logger.warning(f"Autoselected ortholog {p.ensg_id} among {len(ortholog_results)}, but it MAY NOT have the highest percentage identity!")
                        # TODO: this only takes the first ENS id. Implement perc_id (percentage identity) checks!!

                if len(ortholog_results) == 1 and p.gorth_ortholog_status != "definitive": # p.gorth_ortholog_status != "definitive" to prevent duplicates
                    # this is the ideal case -> gOrth returns only one ortholog id = "definitive ortholog"
                    p.gorth_ortholog_exists = True
                    p.gorth_ortholog_status = "definitive"
                    if p.ensg_id is not None: # move existing ensembl id among id synonyms
                        p.id_synonyms.append(p.ensg_id)
                    p.ensg_id = ortholog_results[0]
            logger.info(f"    elapsed: {timer.get_elapsed_formatted('milliseconds')}")
        
        for p in self.products:
            match p.gorth_ortholog_status:
                case "none":
                    num_no_orthologs +=1
                case "indefinitive":
                    num_queried_orthologs_indefinitive += 1
                case "definitive":
                    num_queried_orthologs_definitive +=1
        
        # print some statistics for the user:
        logger.info(f"Finished gOrth batch ortholog query. Results:")
        logger.info(f"  - num definitive orthologs (only 1x ENSG id): {num_queried_orthologs_definitive}")
        logger.info(f"  - num indefinitive orthologs (more than 1 ENSG ids): {num_queried_orthologs_indefinitive}")
        logger.info(f"  - num no orthologs found: {num_no_orthologs}")

    
    def fetch_ortholog_products(
        self,
        target_organism_taxon_number:int = None,
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
        if self.model_settings.ortholog_organisms is None or self.model_settings.ortholog_organisms == []:
            logger.info(f"No orthologous organisms specified. Skipping product ortholog query.")
            return
        
        logger.info("Started fetching ortholog products.")
        self.timer.set_start_time()

        # find target_organism_taxon_number
        if target_organism_taxon_number is None:
            # attempt search in model settings
            if self.model_settings.target_organism.ncbi_id != -1:
                target_organism_taxon_number = self.model_settings.target_organism.ncbi_id
            else:
                # not target organism taxon number specified
                raise Exception(f"No target organism taxon number (id) was specified for the ortholog search!")

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
                        target_organism_taxon_number = target_organism_taxon_number,
                        model_settings=self.model_settings,
                        refetch=refetch,
                        max_connections=max_connections,
                        req_delay=req_delay,
                        semaphore_connections=semaphore_connections,
                    )
                )
            else:
                third_party_db_files = self.model_settings.get_datafile_paths("ALL")
                third_party_db_urls = self.model_settings.get_datafile_urls("ALL")
                human_ortholog_finder = HumanOrthologFinder(
                    goaf=self.goaf,
                    zfin_filepath=third_party_db_files["ortho_mapping_zfin_human"],
                    zfin_download_url=third_party_db_urls["ortho_mapping_zfin_human"],
                    xenbase_filepath=third_party_db_files["ortho_mapping_xenbase_human"],
                    xenbase_download_url=third_party_db_urls["ortho_mapping_xenbase_human"],
                    mgi_filepath=third_party_db_files["ortho_mapping_mgi_human"],
                    mgi_download_url=third_party_db_urls["ortho_mapping_mgi_human"],
                    rgd_filepath=third_party_db_files["ortho_mapping_rgd_human"],
                    rgd_download_url=third_party_db_urls["ortho_mapping_rgd_human"]
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
        target_organism_taxon_number:int = None,
        model_settings:ModelSettings = None,
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
        third_party_db_files = self.model_settings.get_datafile_paths("ALL")
        third_party_db_urls = self.model_settings.get_datafile_urls("ALL")
        human_ortholog_finder = HumanOrthologFinder(
            goaf=self.goaf,
            zfin_filepath=third_party_db_files["ortho_mapping_zfin_human"],
            zfin_download_url=third_party_db_urls["ortho_mapping_zfin_human"],
            xenbase_filepath=third_party_db_files["ortho_mapping_xenbase_human"],
            xenbase_download_url=third_party_db_urls["ortho_mapping_xenbase_human"],
            mgi_filepath=third_party_db_files["ortho_mapping_mgi_human"],
            mgi_download_url=third_party_db_urls["ortho_mapping_mgi_human"],
            rgd_filepath=third_party_db_files["ortho_mapping_rgd_human"],
            rgd_download_url=third_party_db_urls["ortho_mapping_rgd_human"]
        )
        uniprot_api = UniProtApi()
        ensembl_api = EnsemblApi()
        ensembl_api.async_request_sleep_delay = req_delay
        uniprot_api.async_request_sleep_delay = req_delay
        # goaf = GOAnnotationsFile(third_party_db_files['goaf_filepath'])

        # TODO: implement gOrth call -> read orthologs into self.products

        connector = aiohttp.TCPConnector(limit=max_connections, limit_per_host=max_connections)
        semaphore = asyncio.Semaphore(semaphore_connections)
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = []
            for product in self.products:
                if product.had_orthologs_computed is False or refetch is True:
                    # task = product.fetch_ortholog_async(session, human_ortholog_finder, uniprot_api, ensembl_api)
                    task = product.fetch_ortholog_async_semaphore(
                        session,
                        semaphore,
                        self.goaf,
                        target_organism_taxon_number=target_organism_taxon_number,
                        model_settings=model_settings,
                        human_ortholog_finder=human_ortholog_finder,
                        uniprot_api=uniprot_api,
                        ensembl_api=ensembl_api,
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
        """
        Method algorithm:
          - (1) create a dictionary 'reverse_genename_products', where:
                    key = gene name
                    value = list of associated Product objects (with the gene name)
          - (2) iterate 'reverse_genename_products'
                    if more than 1 Product is found in a list of Products mapped to genename,
                    remove all Products from (ReverseLookup).products, create a new Product instance
                    based on the first Product found in the list of 'reverse_genename_products' associated to the current genename.
        """
        logger.info("Started pruning products.")
        start_prod_count = len(self.products)
        self.timer.set_start_time()

        # Create a dictionary that maps genename to a list of products
        reverse_genename_products = {}
        for product in self.products:
            if product.genename is not None:
                # create a mapping to gene name 
                reverse_genename_products.setdefault(product.genename, []).append(product)

        # For each ENSG that has more than one product associated with it, create a new product with all the synonyms
        # and remove the individual products from the list
        for genename, product_list in reverse_genename_products.items():
            if len(product_list) > 1:
                combined_product = Product(id_synonyms=[])
                id_synonyms = []
                for product in product_list:
                    combined_product.update(product)
                    id_synonyms.extend(product.id_synonyms)
                    self.products.remove(product)
                
                # prevent duplicates
                id_synonyms = set(id_synonyms)
                id_synonyms = list(id_synonyms)

                # Create a new product with the collected information and add it to the product list
                self.products.append(combined_product)
                """
                self.products.append(
                    Product(
                        id_synonyms=id_synonyms,
                        taxon=product_list[0].taxon,
                        genename=product_list[0].genename,
                        uniprot_id=product_list[0].uniprot_id,
                        description=product_list[0].description,
                        ensg_id=product_list[0].ensg_id,
                        enst_id=product_list[0].enst_id,
                        refseq_nt_id=product_list[0].refseq_nt_id,
                        mRNA=product_list[0].mRNA,
                        scores=product_list[0].scores,
                        had_orthologs_computed=product_list[0].had_orthologs_computed,
                        had_fetch_info_computed=product_list[0].had_fetch_info_computed,
                        gorth_ortholog_exists=product_list[0].gorth_ortholog_exists
                    )
                )
                """
        
        end_prod_count = len(self.products)
        logger.info(f"Completed product prune operation. Pruned {end_prod_count - start_prod_count} products. Start product count: {start_prod_count} -> End product count: {end_prod_count}")

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
    
    def bulk_ens_to_genename_mapping(self):
        """
        Uses Ensembl's batch mapping using a POST request (https://rest.ensembl.org/documentation/info/lookup_post) to query
        the genename of all the genes with a valid ensembl id.
        """
        # loop over all products, find any ens ids
        ens_ids = []
        ensid_to_product_dict = {}
        for product in self.products:
            if product.ensg_id is not None:
                ens_ids.append(product.ensg_id)
                ensid_to_product_dict[product.ensg_id] = product
            else:
                for id_syn in product.id_synonyms:
                    if "ENS" in id_syn:
                        ens_ids.append(id_syn)
                        ensid_to_product_dict[id_syn] = product
                        break
        
        # fetch
        ensapi = EnsemblApi()
        lookup_result = EnsemblApi.batch_ensembl_lookup(ensapi, ids = ens_ids)
        num_results = len(lookup_result)

        # process results
        ensid_to_genename_dict = {}
        for ens_id, ens_data in lookup_result.items():
            if ens_data is None:
                continue
            genename = ens_data.get('display_name', None)
            if genename is not None:
                ensid_to_genename_dict[ens_id] = genename
        num_genenames = len(ensid_to_genename_dict)
        
        # modify products
        number_updated_genenames = 0 # number of updated genename
        for ens_id, genename in ensid_to_genename_dict.items():
            if ens_id in ensid_to_product_dict:
                p = ensid_to_product_dict[ens_id]
                if p.genename is None:
                    number_updated_genenames += 1
                p.genename = genename
        
        logger.info(f"Finished batch Ensembl to genename mapping. Out of {num_results} mappings, parsed a total of {num_genenames} genenames. Number of updated gene names: {number_updated_genenames}")

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
            for product in tqdm(self.products, "Scoring products"):  # each Product has a field scores - a dictionary between a name of the scoring algorithm and it's corresponding score
                for _score_class in score_classes:
                    # NOTE: Current miRNA scoring (self.score_miRNAs) performs miRNA scoring holistically - in one call for all miRNAs in self.miRNAs. It is pointless to call this function here, as it needs to
                    # be called only once. Here, a function for miRNA scoring has to be called, which displays the top N miRNAs, which bind to the specific product.
                    #
                    # if isinstance(_score_class, basic_mirna_score):
                    #    self.score_miRNAs(_score_class, recalculate=recalculate)
                    #    continue
                    if isinstance(_score_class, basic_mirna_score):
                        continue # just continue, see explanation above

                    if _score_class.name in product.scores and recalculate is True:  # if score already exists and recalculate is set to True
                        product.scores[_score_class.name] = _score_class.metric(product)  # create a dictionary between the scoring algorithm name and it's score for current product
                    elif _score_class.name not in product.scores:  # if score doesn't exist yet
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
                    for SOI in self.target_SOIs:
                        for direction in ["+", "-"]:
                            if "error" in product.scores[_score_class.name][f"{SOI['SOI']}{direction}"]:  # check if there is "error" key
                                continue
                            p_values.append(product.scores[_score_class.name][f"{SOI['SOI']}{direction}"]["pvalue"])
                
                # apply multiple testing correction
                if len(p_values) > 0:
                    reject, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method=self.model_settings.multiple_correction_method)
                    for product in self.products:
                        for SOI in self.target_SOIs:
                            for direction in ["+", "-"]:
                                if "error" in product.scores[_score_class.name][f"{SOI['SOI']}{direction}"]:  # check if there is "error" key
                                    continue
                                product.scores[_score_class.name][f"{SOI['SOI']}{direction}"]["pvalue_corr"] = p_corrected[i]
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

    def get_all_goterms_for_SOI(self, _SOI: str) -> List[GOTerm]:
        """
        Loops through all GO Term objects in self.goterms (initialised from input.txt or from load_model at object creation)
        and adds each GO Term instance to a result list, if any of it's SOIs (goterm.SOIs) are involved in the parameter 'SOI'.

        Returns:
          - List[GOTerm]: a list of all GO Term objects, which are associated with the SOI.

        Example: if SOI = "diabetes", then it will return a list of all the diabetes-associated GO Terms you specified
        in input.txt file, irrespective of direction (either +, - or 0)
        """
        goterms_list = []
        for goterm in self.goterms:
            if any(SOI["SOI"] == _SOI for SOI in goterm.SOIs):
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

    def get_product(self, identifier, identifier_type=None) -> Product:
        """
        Return product based on any id

        Params:
          - identifier: the value to look for
          - identifier_type: specify either 'id_synonyms', 'genename', 'description', 'uniprot_id', 'enst_id', 'refseq_nt_id' or 'mRNA' to speed up search
        """
        ALLOWED_TYPES = ['id_synonyms', 'genename', 'description', 'uniprot_id', 'ensg_id', 'enst_id', 'refseq_nt_id', 'mRNA']
        if identifier_type is not None and identifier_type in ALLOWED_TYPES:
            predictive_dict_name = f"get_product_{identifier_type}_to_product"
            predictive_dict = None
            
            if hasattr(self, predictive_dict_name): # check if predictive dict already exists
                predictive_dict = getattr(self, predictive_dict_name)
            
            if predictive_dict is None: # create predictive dict for this identifier
                res = {}
                for p in self.products:
                    attr = getattr(p, identifier_type)
                    if attr is not None:
                        # exception: add all id synonyms
                        if identifier_type == "id_synonyms": # attr is a list of id synonyms
                            for id_syn in attr:
                                # add both a split and a full id syn
                                if ":" in id_syn:
                                    id_syn_split = id_syn.split(":")[1]
                                    res[id_syn_split] = p
                                res[id_syn] = p  
                            continue
                        # add the rest non-list variables normally
                        res[attr] = p

                setattr(self, predictive_dict_name, res)
                predictive_dict = res
            
            if identifier in predictive_dict:
                return predictive_dict[identifier]
            return None

        for product in self.products:
            if any(identifier in id for id in product.id_synonyms):
                return product
            if any(
                identifier in getattr(product, attr) for attr in [
                    "genename",
                    "description",
                    "uniprot_id",
                    "ensg_id",
                    "enst_id",
                    "refseq_nt_id",
                    "mRNA"
                ] if getattr(product,attr) is not None
            ):
                return product
        logger.debug(f"Couldn't find product for {identifier}")
        return None

    def save_model(self, filepath:str, use_dest_dir:bool=False) -> None:
        """
        Saves the model.
        
        Params:
          - (str) filepath: The relative filepath (ending in .json) to the JSON output file.
          - (bool) use_dest_dir: Whether the ReverseLookup model's destination_dir should be used as the reference for the relative filepath.
                                 Production code should set 'use_dest_dir' to True, so that output files will be saved relatively to the specified destination directory.
                                 Developers should set 'use_dest_dir' to False, so that the output files are saved relatively to the development directory
        """
        self.total_goterms = len(self.goterms)
        self.total_products = len(self.products)

        if ".json" not in filepath:
            filepath = f"{filepath}.json"
            
        if use_dest_dir:
            # use destination dir in project settings - this should be set to True for production-ready code
            # for development, set 'use_dest_dir' to False, so files are saved relatively to filepath branching out from project root directory.
            dest_dir = f"{self.model_settings.destination_dir}"
            filepath = os.path.join(dest_dir, filepath).replace("\\","/")
        
        FileUtil.check_path(filepath)
            
        data = {}
        data["input_filepath"] = self.input_filepath
        data["GO_api_version"] = self.GO_api_version
        data["OBO_version_info"] = self.OBO_version_info
        data["target_SOIs"] = self.target_SOIs
        data["defined_SOIs"] = self.defined_SOIs
        data["go_categories"] = self.go_categories
        data["model_settings"] = self.model_settings.to_json()
        data["miRNA_overlap_treshold"] = self.miRNA_overlap_treshold
        data["execution_times"] = self.execution_times

        # save goterms
        for goterm in self.goterms:
            data.setdefault("goterms", []).append(goterm.__dict__)
        # save products
        for product in self.products:
            data.setdefault("products", []).append(product.__dict__)
        # save stat relevant products
        data["statistically_relevant_products"] = self.statistically_relevant_products
        # save miRNAs
        for mirna in self.miRNAs:
            data.setdefault("miRNAs", []).append(mirna.__dict__)

        logger.info(f"Saving model to: {filepath}")
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

    def perform_statistical_analysis(self, test_name:str="fisher_test", filepath:str="", use_dest_dir:bool = False, two_tailed:bool=False):
        """
        Finds the statistically relevant products, saves them to 'filepath' (if it is provided) and returns a JSON object with the results.

        Parameters:
          - (str) test_name: The name of the statistical test to use for product analysis. It must be either 'fisher_test' (the results of the fisher's test are then used)
          or 'binomial_test' (the results of the binom test are used).
          - (str) filepath: The path to the output file
          - (bool) use_dest_dir: Whether the ReverseLookup model's destination_dir should be used as the reference for the relative filepath.
                                 Production code should set 'use_dest_dir' to True, so that output files will be saved relatively to the specified destination directory.
                                 Developers should set 'use_dest_dir' to False, so that the output files are saved relatively to the development directory
          - (bool) two_tailed: If True, will also analyze the second tail of significance (e.g. first tail is significance of target SOI and insignificance of reverse SOI, 
                               second tail is insignificance of target SOI and significance of reverse SOI)

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
            SOI_PAIR_CODE: [        "diabetes+:angio+": [
                PRODUCT1_DICT               { // product info: id_synonyms, genename, description, ...},
                PRODUCT2_DICT               { // product info: id_synonyms, genename, description, ...},
                ...                         ...
            ],                          ],
            ...                         "diabetes+:obesity+": [...],
                                        "angio+:obesity+": [...]
            }                           }

        The genes (products) for each SOI pair are sorted according to the sum of the p-values, with products with the lowest pvalues (highest
        statistical probabilities) appearing first in the sorted dictionary.

        TODO: implement binomial score, maybe even adv_score and nterms for backwards compatibility
        """

        def sorting_key(product):
            """
            Sorting key used for the sorting of JSON data based on ascending pvalues.
            Fisher test MUST be calculated for this to work.
            """
            pvalue_sum = 0
            for SOI in self.target_SOIs:
                pvalue_SOI = product["scores"][test_name][f"{SOI['SOI']}{SOI['direction']}"]["pvalue_corr"]
                pvalue_sum += pvalue_SOI
            return pvalue_sum
        
        def determine_product_statistical_significance(product:Product, test_name, SOIs, second_tail:bool=False):  
            """
            second_tail: if True, will return True for:
              - if excluding opposite regulation direction check: genes, which statisfy pvalue > (1-p), if len(SOIs) = 1
              - if not excluding opposite regulation direction check: genes, which satisfy pvalue > p (for target SOIs) and pvalue < p (for reverse SOIs)
            
            Returns dictionary:
            {
            'significance': True or False, whether the gene was found to be significant
            'SOIs': the SOIs for which the gene was found to be significant, e.g. chronic_inflammation+:cancer+; obtained by example code: SOI1['SOI']}{SOI1['direction']}:{SOI2['SOI']}{SOI2['direction']
            }
            """
            def get_SOI_groups_label(SOIs, change_direction:bool=False): # example return: diabetes+:angio+
                printSOIs = ""
                for SOI in SOIs:
                    dir = SOI['direction']
                    soi = SOI['SOI']
                    if change_direction:
                        dir = '+' if SOI['direction'] == '-' else '-'
                    printSOIs = f"{printSOIs}:{soi}{dir}"
                s = printSOIs[1: ] # remove first ':'
                return s
            
            # determine if there is only one target SOI defined without complementary reverse SOI,
            # or if there are target SOIs and their complementary reverse SOIs
            has_complementary_SOIs = False
            for SOI in self.target_SOIs:
                soi_name = SOI['SOI']
                direction = 1 if SOI['direction'] == "+" else -1
                for dSOI in self.defined_SOIs:
                    dsoi_name = dSOI['SOI']
                    dsoi_direction = 1 if dSOI['direction'] == "+" else -1
                    if (soi_name == dsoi_name) and (direction + dsoi_direction == 0):
                        has_complementary_SOIs = True
                
            if self.model_settings.exclude_opposite_regulation_direction_check == False and has_complementary_SOIs: # checks both target and the opposite SOIs
                if second_tail: # second tail checks insignificance of target SOI and significance of reverse SOI
                    try:
                        if all( # check target SOI > 0.05
                            float(product.scores[test_name][f"{SOI['SOI']}{SOI['direction']}"]["pvalue_corr"])
                            >= self.model_settings.pvalue
                            for SOI in SOIs   
                        ) and all( # check reverse SOI < 0.05
                            float(product.scores[test_name][f"{SOI['SOI']}{'+' if SOI['direction'] == '-' else '-'}"]["pvalue_corr"])
                            < self.model_settings.pvalue
                            for SOI in SOIs
                        ):
                            return {
                                'significance':True,
                                'SOIs':get_SOI_groups_label(SOIs, change_direction=True)
                            }
                    except KeyError:
                        return {'significance':False}
                # first tail: checks significance of target SOI and insignificance of reverse SOI
                try:
                    if all( # check target SOI < 0.05
                        float(product.scores[test_name][f"{SOI['SOI']}{SOI['direction']}"]["pvalue_corr"])
                        < self.model_settings.pvalue
                        for SOI in SOIs
                    ) and all( # check reverse SOI > 0.05
                        float(product.scores[test_name][f"{SOI['SOI']}{'+' if SOI['direction'] == '-' else '-'}"]["pvalue_corr"])
                        >= self.model_settings.pvalue
                        for SOI in SOIs
                    ):
                        return {
                            'significance':True,
                            'SOIs': get_SOI_groups_label(SOIs, change_direction=False)
                        }
                except KeyError: # most possibly pvalue_corr is missing (can happen if element in contingency table is negative)
                    return {'significance':False}
            else: # checks only target SOI without reverse SOI
                if second_tail: # second tail: check significance of reverse SOI (target SOI >= (1-0.05))
                    try:
                        if all( 
                            float(product.scores[test_name][f"{SOI['SOI']}{SOI['direction']}"].get("pvalue_corr", 0))
                            >= (1-self.model_settings.pvalue)
                            for SOI in SOIs
                        ):
                            return {
                                'significance':True,
                                'SOIs': get_SOI_groups_label(SOIs, change_direction=True)
                            }
                    except KeyError:
                        return {'significance':False}
                try: # first tail: check significance of target SOI (target SOI < 0.05)
                    if all( 
                        float(product.scores[test_name][f"{SOI['SOI']}{SOI['direction']}"].get("pvalue_corr", 1))
                        < self.model_settings.pvalue
                        for SOI in SOIs
                    ):
                        return {
                            'significance':True,
                            'SOIs': get_SOI_groups_label(SOIs, change_direction=False)
                        }
                except KeyError: # most possibly pvalue_corr is missing (can happen if element in contingency table is negative)
                    return {'significance':False}
            return {'significance':False}

        statistically_relevant_products = []  # a list of lists; each member is [product, "SOI1_name_direction:SOI2_name_direction"]
        
        if use_dest_dir:
            # use destination dir in project settings - this should be set to True for production-ready code
            # for development, set 'use_dest_dir' to False, so files are saved relatively to filepath branching out from project root directory.
            dest_dir = f"{self.model_settings.destination_dir}"
            filepath = os.path.join(dest_dir, filepath)
            filepath = filepath.replace("\\", "/")

        for product in self.products:
            # example - given three SOIs: diabetes, angio, obesity, this code iterates through each 2-member combination possible
            #
            # loop iteration \ SOI           diabetes    angio   obesity
            # it. 0  (i=0,j=1)                  |         |
            # it. 1  (i=0,j=2)                  |                  |
            # it. 2  (i=1,j=2)                            |        |
            # j = 3 -> loop condition not met
            # i = 2 -> loop condition not met (i would start on 'obesity', wouldn't find matching pair with j)
            #
            # Each member pair is used to assess statistically relevant genes, which either positively or
            # negatively regulate both of the SOIs in the pair.

            if len(self.target_SOIs) == 1:
                r = determine_product_statistical_significance(
                    product = product,
                    test_name = test_name,
                    SOIs = self.target_SOIs,
                    second_tail=two_tailed
                )
                if r.get('significance'):
                    statistically_relevant_products.append(
                                [
                                    product,
                                    r.get('SOIs')
                                ]
                            )                
                continue # do not advance to the multiple target SOI scoring phase

            # multiple target SOIs scoring phase
            for i in range(len(self.target_SOIs) - 1): # if only 1 target SOI is specified, the whole loop skips because of this condition
                for j in range(i + 1, len(self.target_SOIs)):
                    SOI1 = self.target_SOIs[i]
                    SOI2 = self.target_SOIs[j]
                    pair = [SOI1, SOI2]
                    
                    r = determine_product_statistical_significance(
                        product = product,
                        test_name = test_name,
                        SOIs = pair,
                        second_tail=two_tailed
                    )
                    if r.get('significance'):
                        statistically_relevant_products.append(
                                [
                                    product,
                                    r.get('SOIs')
                                ]
                            )
        
        statistically_relevant_products_final = {} # dictionary between two SOIs (eg. angio+:diabetes+) and all statistically relevant products or a single target SOI (eg. angio+) and all statistically relevant products
        if len(self.target_SOIs) == 1:
            # do not advance to the multiple target SOI score analysis phase
            for element in statistically_relevant_products:
                prod = element[0]
                SOI_label = element[1]
                if SOI_label not in statistically_relevant_products_final:
                    statistically_relevant_products_final[SOI_label] = []
                statistically_relevant_products_final[SOI_label].append(prod.__dict__)
            
            logger.info(f"Displaying significant genes:")
            statistically_relevant_products_final_sorted = {}
            num_significant_genes = 0
            for SOIs_label,significant_products in statistically_relevant_products_final.items():
                statistically_relevant_products_for_SOIs_label_sorted = sorted( # sorts each subgroup (each SOI label)
                    significant_products,
                    key= lambda gene: sorting_key(gene)
                    )
                statistically_relevant_products_final_sorted[SOIs_label] = statistically_relevant_products_for_SOIs_label_sorted
                # update gene counts
                num_significant_genes += len(statistically_relevant_products_for_SOIs_label_sorted)
                logger.info(f"  - {SOIs_label} : {len(statistically_relevant_products_for_SOIs_label_sorted)} genes")
            logger.info(f"  - total significant genes: {num_significant_genes}")
            
            # update significant genes
            self.statistically_relevant_products = statistically_relevant_products_final_sorted            
            
            logger.info(f"Saving statistically relevant products to {os.path.abspath(filepath)}")
            JsonUtil.save_json(
                data_dictionary=statistically_relevant_products_final_sorted,
                filepath=filepath
            )
            return statistically_relevant_products_final_sorted

        # * multiple target SOIs score analysis phase *
        # statistically_relevant_products stores a list of lists, each member list is a Product object bound to a specific pair code (e.g. angio+:diabetes+).
        # statistically_relevant_products_final is a dictionary. It's keys are SOIs pair codes (e.g. angio+:diabetes+), each key holds a list of all statistically relevant products for the SOI pair
        # (eg. if angio+:diabetes+ it holds all products, which positively regulate both angiogenesis and diabetes)
        SOI_pairs = [] # each element is a code binding two SOIs and their direction, eg. angio+:diabetes+
        for i in range(len(self.target_SOIs) - 1):
            for j in range(i + 1, len(self.target_SOIs)):
                SOI1 = self.target_SOIs[i]
                SOI2 = self.target_SOIs[j]
                pair_code = f"{SOI1['SOI']}{SOI1['direction']}:{SOI2['SOI']}{SOI2['direction']}"
                SOI_pairs.append(pair_code)
                statistically_relevant_products_final[pair_code] = []  # initialise to empty list

        for element in statistically_relevant_products:
            # each element is a list [product, "SOI1_name_direction:SOI2_name_direction"]
            prod = element[0]
            SOI_pair_code = element[1]
            if SOI_pair_code not in statistically_relevant_products_final: # takes care also of the SOI pair codes if two-tailed test was used (this changes the SOI directions)
                statistically_relevant_products_final[SOI_pair_code] = []
            statistically_relevant_products_final[SOI_pair_code].append(prod.__dict__)
        
        logger.info(f"Displaying significant genes:")
        num_significant_genes = 0
        # sort the genes based on the ascending sum of pvalues (lowest pvalues first)
        statistically_relevant_products_final_sorted = {}
        for SOIs_label,significant_products in statistically_relevant_products_final.items():
            statistically_relevant_products_for_SOIs_label_sorted = sorted(
                significant_products,
                key= lambda gene: sorting_key(gene)
            )
            statistically_relevant_products_final_sorted[SOIs_label] = statistically_relevant_products_for_SOIs_label_sorted
            # update gene counts
            num_significant_genes += len(statistically_relevant_products_for_SOIs_label_sorted)
            logger.info(f"  - {SOIs_label} : {len(statistically_relevant_products_for_SOIs_label_sorted)} genes")
        logger.info(f"  - total significant genes: {num_significant_genes}")

        """
        for i in range(len(self.target_SOIs) - 1):
            for j in range(i + 1, len(self.target_SOIs)):
                SOI1 = self.target_SOIs[i]
                SOI2 = self.target_SOIs[j]
                pair_code = f"{SOI1['SOI']}{SOI1['direction']}:{SOI2['SOI']}{SOI2['direction']}"
                statistically_relevant_products_for_SOI_pair = statistically_relevant_products_final[pair_code]
                statistically_relevant_products_for_SOI_pair_sorted = sorted(
                    statistically_relevant_products_for_SOI_pair,
                    key=lambda gene: sorting_key(gene)
                )
                statistically_relevant_products_final_sorted[pair_code] = statistically_relevant_products_for_SOI_pair_sorted
        """
        
        # TODO: save statistical analysis as a part of the model's json and load it up on startup
        self.statistically_relevant_products = statistically_relevant_products_final_sorted
        logger.info(f"Finished with product statistical analysis. Found {len(statistically_relevant_products)} statistically relevant products. p = {self.model_settings.pvalue}")
        JsonUtil.save_json(
            data_dictionary=statistically_relevant_products_final_sorted,
            filepath=filepath
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
    def load_model(cls, filepath: str, destination_dir:str = None) -> "ReverseLookup":
        """
        Loads the model representation from an existing .json file.
        
        Params:
          - (str) filepath: the path to the .json file of a previously saved model
          - (str) destination_dir: The destination_dir should be used in production code, especially if the ReverseLookup model is also being saved using a use_dest_dir = True parameter.
        
        This method also parses GO Term parents for GO Terms, if setting include_indirect_annotations is True.
        """        
        if destination_dir is not None:
            # use destination dir in project settings - this should be set to True for production-ready code
            # for development, set 'use_dest_dir' to False, so files are saved relatively to filepath branching out from project root directory.
            filepath = os.path.join(destination_dir, filepath)
            filepath = filepath.replace("\\", "/")
        
        logger.info(f"Loading ReverseLookup model from {os.path.abspath(filepath)}")
        data = JsonUtil.load_json(filepath)
        if data == {}:
            logger.warning(f"Data is EMPTY!")
            
        target_SOIs = data["target_SOIs"]
        defined_SOIs = data["defined_SOIs"]
        miRNA_overlap_treshold = data["miRNA_overlap_treshold"]
        
        input_filepath = data.get("input_filepath", None)
        GO_api_version = data.get("GO_api_version", None)
        OBO_version_info = data.get("OBO_version_info", None)

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
        if destination_dir is not None:
            settings.destination_dir = destination_dir

        goterms = []
        for goterm_dict in data["goterms"]:
            goterms.append(GOTerm.from_dict(goterm_dict))

        obo_parser = None
        if settings.include_indirect_annotations is True:
            if settings.datafile_paths != {} and "go_obo" in settings.datafile_paths:
                if(
                    settings.datafile_paths['go_obo'] is not None
                    and settings.get_datafile_path("go_obo") != ""
                ):
                    obo_parser = OboParser(obo_filepath=settings.get_datafile_path("go_obo"), obo_download_url=settings.get_datafile_url("go_obo"))
            else:
                obo_parser = OboParser()
            for goterm in goterms:
                assert isinstance(goterm, GOTerm) # TODO: FIX HERE !!!!!!!!! obo now returns a dict and not a goterm!!!
                if goterm.parent_term_ids == [] or goterm.parent_term_ids is None:
                    goterm_obo = GOTerm.from_dict(obo_parser.all_goterms[goterm.id].__dict__)  # obo representation of this goterm is in json form
                    goterm.update(goterm_obo)  # update current goterm with information from .obo file

                    goterm_parent_ids = obo_parser.get_parent_terms(goterm.id)  # calculate parent term ids for this goterm
                    goterm.parent_term_ids = goterm_parent_ids  # update parent term ids

        products = []
        for product_dict in data.get("products", []):
            products.append(Product.from_dict(product_dict))

        miRNAs = []
        for miRNAs_dict in data.get("miRNAs", []):
            miRNAs.append(miRNA.from_dict(miRNAs_dict))
            
        logger.info(f"Loaded ReverseLookup model from {os.path.abspath(filepath)}")
        logger.info(f"  - num goterms: {len(goterms)}")
        logger.info(f"  - num products: {len(products)}")
        logger.info(f"  - num miRNAs: {len(miRNAs)}")
        logger.info(f"  - model_settings: {settings.to_json()}")
        logger.info(f"  - input_filepath: {input_filepath}")
        logger.info(f"  - example goterm (0): {goterms[0].to_json()}")

        return cls(
            goterms,
            target_SOIs,
            products,
            miRNAs,
            miRNA_overlap_treshold,
            execution_times=execution_times,
            statistically_relevant_products=statistically_relevant_products,
            go_categories=go_categories,
            model_settings=settings,
            obo_parser=obo_parser,
            input_filepath=input_filepath,
            GO_api_version=GO_api_version,
            OBO_version_info=OBO_version_info,
            defined_SOIs=defined_SOIs
        )

    @classmethod
    def from_input_file(cls, filepath: str, destination_dir:str=None) -> "ReverseLookup":
        """
        Creates a ReverseLookup object from a text file.

        Args:
            filepath (str): The path to the input text file.
            destination_dir (str): The destination directory where output files are stored

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        # Define constants used in parsing the file
        LINE_ELEMENT_DELIMITER = "\t"  # Data is tab separated
        COMMENT_DELIMITER = "#"  # Character used to denote a comment
        LOGIC_LINE_DELIMITER = "###" # Special set of characters to denote a "logic line"

        target_SOIs = []
        go_categories = []
        go_terms = []
        settings = ModelSettings()
        settings.destination_dir = destination_dir
        
        def process_comment(line):
            """
            Processes a comment in the line: returns the part of the line before the comment. The input file should be structured to contain
            three sections - 'settings', 'states_of_interest' and 'GO_terms', annotated using the LOGIC_LINE_DELIMITER.

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
        
        def process_evidence_codes(evidence_code_instructions:str, all_evidence_codes:dict):
            """
            Processes an evidence code instruction, compares the evidence codes to 'all_evidence_codes' and determines which evidence codes should be deemed as valid for the research.
            Note that all_evidence_codes MUST be already computed at the point of calling this function. The result of this function is the creation of
            a list of valid evidence codes.

            Example evidence_code_instructions: "author_statement(~),curator_statement(IC)"

            For more information, refer to the explanation in the demo input.txt file
            """
            if all_evidence_codes is None or all_evidence_codes == {}:
                raise Exception("all_evidence_codes in process_evidence_codes function is None or {}. Make sure to supply evidence_code section in the input.txt BEFORE the settings section.")
            
            # create a list of evidence code groups instructions (multiple or a single)
            evidence_code_instructions = evidence_code_instructions.split(",") if "," in evidence_code_instructions else [evidence_code_instructions]
            # temporary dictionary to store processing results - so the program knows which groups of evidence codes have been processed
            # if any group hasn't been processed (meaning it wasn't set by the user), automatically add all items of the unprocessed group to valid evidence codes
            evidence_code_groups_processing_states = {}
            for code_group,evidence_codes in all_evidence_codes.items():
                evidence_code_groups_processing_states[code_group] = False
            
            # process instruction evidence code
            valid_evidence_codes = []
            for instruction in evidence_code_instructions: # example instruction: author_statement(~)
                instruction_code_group = instruction.split("(")[0]
                instruction_evidence_codes = instruction.split("(")[1].split(")")[0]
                negate_group = False # if '!' is used to negate the meaning (exclude the group or evidence codes)
                if '!' in instruction_code_group:
                    negate_group = True
                    instruction_code_group = instruction_code_group.replace("!", "")

                evidence_code_groups_processing_states[instruction_code_group] = True

                if instruction_evidence_codes == "~":
                    if negate_group == False:
                        # valid_evidence_codes += all_evidence_codes[instruction_code_group]
                        for full_code in all_evidence_codes[instruction_code_group]:
                            eco_evidence_code_id = full_code.split("_")[1]
                            valid_evidence_codes.append(eco_evidence_code_id)
                else: # instruction_evidence_codes contains specific evidence codes
                    instruction_evidence_codes = instruction_evidence_codes.split(",") if "," in instruction_evidence_codes else [instruction_evidence_codes]
                    if negate_group == False:
                        # experimental(EXP) = include only EXP from experimental codes, exclude the rest of the experimental codes
                        for instruction_evidence_code in instruction_evidence_codes:
                            # instruction_evidence_code is only EXP or IBA (not the ECO:xxxx identifier) -> convert to ECO identifier
                            valid_evidence_codes.append(settings.evidence_codes_to_ecoids.get(instruction_evidence_code))
                    else: # !experimental(EXP) = exclude only EXP from experimental codes, include the rest
                        for full_evidence_code in all_evidence_codes[instruction_code_group]:
                            for instruction_evidence_code in instruction_evidence_codes:
                                if instruction_evidence_code not in full_evidence_code and full_evidence_code.split("_")[1] not in valid_evidence_codes:
                                    valid_evidence_codes.append(full_evidence_code.split("_")[1])
            
            # add all non-user-specified groups
            for code_group,val in evidence_code_groups_processing_states.items():
                if val == False: # if code group wasn't processed
                    for full_code in all_evidence_codes[code_group]:
                        eco_evidence_code_id = full_code.split("_")[1]
                        valid_evidence_codes.append(eco_evidence_code_id)
                    #valid_evidence_codes += all_evidence_codes[code_group]
            
            return valid_evidence_codes

        filepath_readlines = 0
        with open(filepath, 'r') as f:
            filepath_readlines = len(f.readlines())

        def process_file(filepath: str):
            all_evidence_codes = {}
            evidence_codes_to_ecoids = {} # maps evidence codes (e.g. EXP) to respective ECO ids (e.g. ECO:0000269)
            with open(filepath, "r") as read_content:
                # read_lines = read_content.read().splitlines()[2:]  # skip first 2 lines
                read_lines = read_content.read().splitlines()
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
                    elif "states_of_interest" in line:
                        section = "states_of_interest"
                        continue
                    elif "categories" in line:
                        section = "categories"
                        continue
                    elif "GO_terms" in line:
                        section = "GO"
                        continue
                    if "evidence_code_groups" in line:
                        section = "evidence_code_groups"
                        continue
                    
                    if section == "settings":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        if len(chunks) < 2:
                            # means there is no setting value, e.g. only "ortholog_organisms" is specified without any ortholog organisms specified
                            continue
                        setting_name = chunks[0]
                        setting_value = chunks[1]  # is string now
                        setting_optionals = None
                        if len(chunks) > 2: # this means that an optional setting value was specified eg. 'include_indirect_annotations   True    c'
                            setting_optionals = chunks[2]
                            
                        if setting_value == "True" or setting_value == "true":
                            setting_value = True
                        if setting_value == "False" or setting_value == "false":
                            setting_value = False
                        if setting_name == "pvalue":
                            setting_value = float(chunks[1])
                        if setting_name == "multiple_correction_method":
                            logger.info(f"Using {setting_value} as multiple correction.")

                        if setting_name == "include_indirect_annotations" and setting_optionals is not None:
                            # this means that a setting optional of either parent or children direction of indirect annotations was specified
                            # eg. 'include_indirect_annotations   True    c'
                            if len(setting_optionals) == 1:
                                settings.indirect_annotations_direction = setting_optionals
                            else:
                                raise Exception(f"Setting optionals value for include_indirect_annotations is not specified correctly. It must be either 'c' or 'p'.")
                                
                        if setting_name == "goterms_set":
                            if setting_value != 'all':
                                if ',' in setting_value:
                                    # e.g. split 'human,rattus_norvegicus' into ["human", "rattus_norvegicus"]
                                    setting_value = setting_value.split(',')
                                else:
                                    # 'human' -> ["human"]
                                    setting_value = [setting_value]

                        if setting_name == "target_organism":
                            organism_info = OrganismInfo.parse_organism_info_str(metadata=setting_value)
                            setting_value = organism_info

                        if setting_name == "ortholog_organisms":
                            organism_info_dict = {}
                            for organism_info_str in setting_value.split(","): # split at commas
                                organism_info = OrganismInfo.parse_organism_info_str(metadata=organism_info_str)
                                # create multiple annotations in dict both for the label and for the ncbitaxon full id
                                if organism_info.label != "":
                                    organism_info_dict[organism_info.label] = organism_info
                                if organism_info.ncbi_id_full != "":
                                    organism_info_dict[organism_info.ncbi_id_full] = organism_info
                            setting_value = organism_info_dict
                        
                        # finally, set the setting
                        settings.set_setting(setting_name=setting_name, setting_value=setting_value)

                        if setting_name == "evidence_codes":
                            # example line: evidence_codes \t experimental(~),phylogenetic(~),computational_analysis(~),author_statement(TAS),curator_statement(IC),!electronic(~)
                            valid_evidence_codes = process_evidence_codes(setting_value, settings.all_evidence_codes)
                            settings.valid_evidence_codes = valid_evidence_codes

                            # display evidence codes in user friendly format
                            valid_true_evidence_codes = []
                            ecoids_to_evidence_codes = DictUtil.reverse_dict(settings.evidence_codes_to_ecoids)
                            for eco_id_code in settings.valid_evidence_codes:
                                valid_true_evidence_codes.append(ecoids_to_evidence_codes.get(eco_id_code))
                            logger.info(f"Using the following evidence codes: {valid_true_evidence_codes}")
                            
                    elif section == "evidence_code_groups":
                        # line example: experimental \t EXP,IDA,IPI,IMP,IGI,IEP,HTP,HDA,HMP,HGI,HEP
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        evidence_code_group = chunks[0]
                        if "," in chunks[1]:
                            evidence_codes = chunks[1].split(",")
                        else: # no comma, thus a single value -> convert to list
                            evidence_codes = [chunks[1]]
                        all_evidence_codes[evidence_code_group] = evidence_codes
                        settings.all_evidence_codes = all_evidence_codes

                        # process eco id
                        for full_code in evidence_codes:
                            if "_" not in full_code:
                                raise Exception(f"The character '_' isn't present in evidence code {full_code}. Is the full evidence code specified in the format CODE_ECOid? For example: EXP_ECO:0000269")
                            true_evidence_code = full_code.split("_")[0] # e.g. EXP
                            eco_id = full_code.split("_")[1] # e.g. ECO:0000269
                            evidence_codes_to_ecoids[true_evidence_code] = eco_id
                        settings.evidence_codes_to_ecoids = evidence_codes_to_ecoids

                    elif section == "filepaths":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        datafile_name = chunks[0]
                        datafile_local_path = chunks[1]
                        if len(chunks) >= 2: # for backwards compatibility
                            datafile_download_url = chunks[2]
                            try:
                                organism = chunks[3]
                            except Exception as e:
                                logger.info(f"'organism' wasn't defined for {datafile_name}. Was that intended?")
                        else:
                            datafile_download_url = None
                            organism = None
                        
                        settings.datafile_paths[datafile_name] = {
                            'organism': organism,
                            'local_filepath': datafile_local_path,
                            'download_url': datafile_download_url
                        }
                    
                    elif section == "states_of_interest":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        target_SOIs.append({"SOI": chunks[0], "direction": chunks[1]})
                    
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
                                "SOIs": {
                                    "SOI": chunks[1],
                                    "direction": chunks[2],
                                },
                                "weight": int(chunks[3]),
                                "name": chunks[4],
                            }
                        else:
                            d = {
                                "id": chunks[0],
                                "SOIs": {
                                    "SOI": chunks[1],
                                    "direction": chunks[2],
                                },
                                "weight": int(chunks[3]),
                            }
                        if not any(d["id"] == goterm.id for goterm in go_terms):  # TODO: check this !!!!!
                            go_terms.append(GOTerm.from_dict(d))
                        else:  # TODO: check this !!!!!
                            next(goterm for goterm in go_terms if d["id"] == goterm.id).add_SOI(
                                {"SOI": chunks[1], "direction": chunks[2]}
                            )

        try:
            process_file(filepath)
        except OSError:
            logger.error(f"ERROR while processing input file at filepath {filepath}")
            return
        
        obo_parser = None
        if settings.include_indirect_annotations:
            if settings.datafile_paths != {} and "go_obo" in settings.datafile_paths:
                if(
                    settings.datafile_paths['go_obo'] is not None
                    and settings.get_datafile_path("go_obo") != ""
                ):
                    obo_parser = OboParser(obo_filepath=settings.get_datafile_path("go_obo"), obo_download_url=settings.get_datafile_url("go_obo"))
            else:
                obo_parser = OboParser()
            logger.info(f"Starting OboParser to find all GO Term parents and children using data file {obo_parser.filepath}")

            # update goterms to include all parents and children
            with logging_redirect_tqdm():
                for goterm in tqdm(go_terms, desc="Compute indirect nodes"):
                    assert isinstance(goterm, GOTerm)
                    if goterm.parent_term_ids == [] or goterm.parent_term_ids is None:
                        goterm_obo = GOTerm.from_dict(obo_parser.all_goterms[goterm.id].__dict__)  # obo representation of this goterm
                        goterm.update(goterm_obo)  # update current goterm with information from .obo file

                        goterm_parent_ids = obo_parser.get_parent_terms(goterm.id)  # calculate parent term ids for this goterm
                        goterm_children_ids = obo_parser.get_child_terms(goterm.id)  # calculdate child term ids for this goterm
                        goterm.parent_term_ids = goterm_parent_ids  # update parent term ids
                        # goterm.child_term_ids = goterm_children_ids  # update child term ids
            logger.info("Indirect annotations have been computed.")
        
        logger.info(f"Computing defined SOIs:")
        defined_SOIs_fullnames = {} # dict linking a SOI full name (e.g. 'chronic_inflammation+') to its dict representation (e.g. {'SOI': "chronic_inflammation", 'direction': "+"})
        for goterm in go_terms:
            assert isinstance(goterm, GOTerm)
            goterm_SOIs = goterm.SOIs
            for goterm_SOI in goterm_SOIs:
                fullname = f"{goterm_SOI['SOI']}{goterm_SOI['direction']}"
                if fullname not in defined_SOIs_fullnames:
                    defined_SOIs_fullnames[fullname] = goterm_SOI
        
        defined_SOIs = []
        for fullname, SOI_dict in defined_SOIs_fullnames.items():
            logger.info(f"  - {fullname}")
            defined_SOIs.append(SOI_dict)
            
        logger.info("Creating model from input file with:")
        logger.info(f"  - input file filepath: {filepath}")
        logger.info(f"  - destination dir: {destination_dir}")
        logger.info(f"  - input file line count: {filepath_readlines}")
        logger.info(f"  - count GO Terms: {len(go_terms)} ")
        logger.info(f"  - target_SOIs: {target_SOIs}")
        logger.info(f"  - GO categories: {go_categories}")
        logger.info(f"  - model settings: {settings.to_json()}")
        logger.info(f"  - obo_parser: {obo_parser}")
        return cls(
            go_terms,
            target_SOIs=target_SOIs,
            defined_SOIs=defined_SOIs,
            go_categories=go_categories,
            model_settings=settings,
            obo_parser=obo_parser,
            input_filepath=filepath
        )

    @classmethod
    def from_dict(cls, data: Dict[str, List[Dict]]) -> "ReverseLookup":
        """
        Creates a ReverseLookup object from a dictionary.

        Args:
            data (dict): A dictionary containing a representation of a ReverseLookup instance

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        goterms = [GOTerm.from_dict(d) for d in data["goterms"]]
        target_SOIs = data["target_SOIs"]
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
        
        input_filepath = data.get("input_filepath", None)

        logger.info("Model creation from dict complete.")
        return cls(
            goterms,
            target_SOIs,
            go_categories=go_categories,
            model_settings=settings,
            input_filepath=input_filepath
        )

    def _debug_shorten_GO_terms(self, count):
        """
        Shortens the amount of GO terms to the specified 'count', for debugging purposes.
        """
        if count < len(self.goterms):
            self.goterms = self.goterms[0:count]

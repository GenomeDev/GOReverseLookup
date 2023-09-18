import aiohttp
import asyncio
from typing import List, Dict, Optional

from ..web_apis.GOApi import GOApi
from ..parse.GOAFParser import GOAnnotationsFile
from ..util.CacheUtil import Cacher

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class GOTerm:
    def __init__(
        self,
        id: str,
        processes: List[Dict] = None,
        name: Optional[str] = None,
        description: Optional[str] = None,
        category: Optional[str] = None,
        parent_term_ids: Optional[List[str]] = None,
        is_obsolete: bool = False,
        weight: float = 1.0,
        products: List[str] = [],
        http_error_codes: dict = {},
    ):
        """
        A class representing a Gene Ontology term.

        Args:
            id (str): The ID of the GO term.
            processes (List[Dict]): [{"process" : "angio", "direction" : "+"}]
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            weight (float): The weight of the GO term.
            products (list): Products associated with the term (optional).
            category (str): biological_process, molecular_activity or cellular_component
            parent_term_ids (list[str]): GO ids of the parent terms (parsed from .obo)
            is_obsolete (bool): if the term is labelled as obsolete in the .obo file
        """
        self.id = id
        self.processes = processes if isinstance(processes, list) else [processes]
        self.name = name
        self.description = description
        if weight is None:  # bugfix for any mistakes
            weight = 1.0
        self.weight = float(weight)
        self.products = products
        self.http_error_codes = {}  # used for errors happening during server querying; for example, a dict pair 'products': "HTTP Error ..." signifies an http error when querying for GO Term's products
        self.category = category
        self.parent_term_ids = parent_term_ids
        self.is_obsolete = is_obsolete

    def update(self, goterm, overwrite_existing: bool = False):
        """
        Updates the member values of this GOTerm instance with the member values of the parameter 'goterm' GOTerm instance.
        This function will automatically check that the ids of this GOTerm and the parameter 'goterm' match before performing any functionality.
        If the ids don't match, en Error will be raised.

        Params:
          - goterm: should be a GOTerm instance
          - overwrite_existing: if True, will overwrite EVERY value of the current GOTerm instance with the parameter 'goterm' values.
                                if False, will only overwrite values of the current GOTerm instance which are None or "" (empty strings)
        """

        def overwrite_attributes(goterm):
            assert isinstance(goterm, GOTerm)
            for attr_name in dir(
                goterm
            ):  # loop through all class attributes of input goterm
                if not callable(
                    getattr(goterm, attr_name)
                ) and not attr_name.startswith("__"):
                    if hasattr(
                        self, attr_name
                    ):  # check if the same attribute exists in self
                        if overwrite_existing is True:
                            setattr(
                                self, getattr(goterm, attr_name)
                            )  # update self attribute with goterm's attribute
                        else:  # overwrite_existing == False
                            attr_value_self = getattr(self, attr_name)
                            if attr_value_self is None or attr_value_self is False:
                                setattr(
                                    self, attr_name, getattr(goterm, attr_name)
                                )  # update self attribute with goterm's attribute only if self attribute is None or False

        assert isinstance(goterm, GOTerm)
        if self.id == goterm.id:  # perform update only if ids match
            overwrite_attributes(goterm=goterm)
        else:
            raise Exception(
                f"Input GOTerm id {goterm.id} doesn't match with current instance"
                f" GOTerm id {self.id}"
            )

    def copy(self):
        """
        Creates a copy of self.

        Example problem: class ReverseLookup holds a list of GOTerms and you want to copy these terms to a new list without storing references to the old terms.
        If you just loop over existing GOTerm instances in ReverseLookup and append them to a new list, changes in one or the other list will affect GOTerms in both lists,
        since you are storing a reference to that object in the list, not creating a copy of the object.

        Solution: use the copy method to create copies of the GOTerm.
            goterm_copies = []
            for goterm in (ReverseLookup).goterms:
                goterm_copies.append(goterm.copy())
        """
        return GOTerm(
            self.id,
            self.processes,
            self.name,
            self.description,
            self.weight,
            self.products,
            self.http_error_codes,
        )

    def fetch_name_description(self, api: GOApi):
        """
        Sets the "name" and "description" member field of the GO Term. The link used to query for the response is http://api.geneontology.org/api/ontology/term/{term_id}.
        This function sets the "name" field of the GO Term to response['label'] and the "description" field to response['definition']

        Parameters:
          - (GOApi) api: a GOApi instance

        Usage and calling:
            api = GOApi()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_name_description(api)

        Example api response is:
        {
            'goid': 'GO:1903589',
            'label': 'positive regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis',
            'definition': 'Any process that activates or increases the frequency, rate or extent of blood vessel endothelial cell proliferation involved in sprouting angiogenesis.',
            'creation_date': '2014-11-04T11:39:47Z',
            'synonyms': [
                'up regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis',
                'up-regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis',
                'upregulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis'
                ],
            'relatedSynonyms': [
                'activation of blood vessel endothelial cell proliferation during sprouting angiogenesis',
                'positive regulation of blood vessel endothelial cell proliferation during sprouting angiogenesis',
                'up regulation of blood vessel endothelial cell proliferation during sprouting angiogenesis',
                'up-regulation of blood vessel endothelial cell proliferation during sprouting angiogenesis',
                'upregulation of blood vessel endothelial cell proliferation during sprouting angiogenesis'
                ],
            'alternativeIds': [''],
            'xrefs': [''],
            'subsets': ['']
        }
        """
        logger.info("Fetching GO Term names (labels) and descriptions (definitions).")
        data = api.get_data(self.id)
        if data:
            if "label" in data:
                self.name = data["label"]
            if "definition" in data:
                self.description = data["definition"]
            logger.info(f"Fetched name and description for GO term {self.id}")

    async def fetch_name_description_async(self, api: GOApi, req_delay=0.1):
        async with aiohttp.ClientSession() as session:
            url = api.get_data(self.id, get_url_only=True)
            response = await session.get(url)
            if response.status == 200:
                data = await response.json()
                await asyncio.sleep(req_delay)
                if "label" in data:
                    self.name = data["label"]
                if "definition" in data:
                    self.description = data["definition"]
                # logger.info(f"Fetched name and description for GO term {self.id}")
                # print out only 15 desc chars not to clutter console
                logger.info(
                    f"GOid {self.id}: name = {self.name}, description ="
                    f" {self.description[:15]}..."
                )
            else:
                logger.info(
                    f"Query for url {url} failed with response code {response.status}"
                )

    def fetch_products(self, source):
        """
        Fetches UniProtKB products associated with a GO Term and sets the "products" member field of the GO Term to a list of all associated products.
        The link used to query for the response is http://api.geneontology.org/api/bioentity/function/{term_id}/genes.

        The product IDs can be of any of the following databases: UniProt, ZFIN, Xenbase, MGI, RGD
        [TODO: enable the user to specify databases himself]

        Parameters:
          - source: can either be a GOApi instance (web-based download) or a GOAnnotationFile isntance (file-based download)

        Usage and calling:
            source = GOApi() OR source = GOAnnotationFile()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_products(source)
        """
        if isinstance(source, GOApi):
            products = source.get_products(self.id)
            if products:
                self.products = products
        elif isinstance(source, GOAnnotationsFile):
            products = source.get_products(self.id)
            if products:
                self.products = products

    async def fetch_products_async_v1(self, api: GOApi, delay=0.0):
        """
        Asynchronously queries the products using api.get_products_async, using the
        await keyword to wait for the request to finish inside the api.

        delay is in seconds
        """
        await asyncio.sleep(delay)
        # products = await api.get_products_async(self.id)
        products = await api.get_products_async_notimeout(
            self.id
        )  # testing variant of the function, comment this later
        if isinstance(products, list) and products != []:
            self.products = products
        elif isinstance(products, str):
            # str reports an HTTP Error code
            if "HTTP Error" in products:
                error_report = products
                self.http_error_codes["products"] = error_report

    async def fetch_products_async_v3(
        self,
        session: aiohttp.ClientSession,
        request_params={"rows": 50000},
        req_delay=0.5,
        max_retries=3,
    ):
        """
        A better variant of get_products_async. Doesn't include timeout in the url request, no retries.
        Doesn't create own ClientSession, but relies on external ClientSession, hence doesn't overload the server as does the get_products_async_notimeout function.

        # Previous algorithm created one aiohttp.ClientSession FOR EACH GOTERM. Therefore, each ClientSession only had one connection,
        # and the checks for connection limiting weren't enforeced. During runtime, there could be as many as 200 (as many as there are goterms)
        # active ClientSessions, each with only one request. You should code in the following manner:
        #
        # async def make_requests():
        #    connector = aiohttp.TCPConnector(limit=20, limit_per_host=20)
        #    async with aiohttp.ClientSession(connector=connector) as session:
        #        urls = [...]  # List of URLs to request
        #        for url in urls:
        #            await asyncio.sleep(1)  # Introduce a 1-second delay between requests
        #            response = await session.get(url)
        #            # Process the response
        """
        # data key is in the format [class_name][function_name][function_params]
        data_key = f"[{self.__class__.__name__}][{self.fetch_products_async_v3.__name__}][go_id={self.id}]"
        previous_data = Cacher.get_data("go", data_key)
        if previous_data is not None and previous_data != []:
            logger.debug(f"Cached previous product fetch data for {self.id}. Data = {previous_data}")
            self.products = previous_data
            return previous_data

        APPROVED_DATABASES = [
            ["UniProtKB", ["NCBITaxon:9606"]],
            ["ZFIN", ["NCBITaxon:7955"]],
            # ["RNAcentral", ["NCBITaxon:9606"]],
            ["Xenbase", ["NCBITaxon:8364"]],
            ["MGI", ["NCBITaxon:10090"]],
            ["RGD", ["NCBITaxon:10116"]],
        ]
        url = f"http://api.geneontology.org/api/bioentity/function/{self.id}/genes"
        params = request_params  # 10k rows resulted in 56 mismatches for querying products for 200 goterms (compared to reference model, loaded from synchronous query data)

        retries = 0
        data = None
        response_url = ""
        for _ in range(max_retries):
            possible_http_error_text = ""
            if retries == (max_retries - 1):
                logger.warning(f"Exceeded max retries when parsing {self.id}")
                raise Exception(
                    f"Exceeded max retries when parsing {self.id}. HTTP error text ="
                    f" {possible_http_error_text}"
                )
            retries += 1

            previous_response = Cacher.get_data("url", url)
            if previous_response is not None:
                data = previous_response
                break  # previous response json was cached, break the loop
            else:
                try:
                    await asyncio.sleep(req_delay)
                    response = await session.get(url, params=params)
                    if (
                        response.status != 200
                    ):  # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                        possible_http_error_text = (
                            f"HTTP Error when parsing {self.id}. Response status ="
                            f" {response.status}"
                        )
                        logger.warning(possible_http_error_text)
                        continue
                    response_url = response.url
                    data = await response.json()
                    if data is not None:
                        Cacher.store_data("url", url, data)
                        logger.debug(f"Cached async product fetch data for {self.id}")
                        break  # reponse json was obtained, break the loop
                except (aiohttp.ClientConnectionError, aiohttp.ClientPayloadError) as e:
                    logger.warning(f"Error when fetching products for {self.id}")
                    logger.warning(f"  - attempted url: {url}")
                    # raise e # TODO: implement logic for errors when querying GO Term products - refetch?

        products_set = set()
        for assoc in data["associations"]:
            if assoc["object"]["id"] == self.id and any(
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
        if products == []:
            logger.warning(
                f"Found no products for GO Term {self.id} (name = {self.name})!"
            )
            logger.warning(f"  - response.url: {response_url}")
            logger.warning(f"  - response.json(): {data}")
            return []
            # if len(data) < 500:
            #    logger.debug(f"Response json: {data}")

        self.products = products
        logger.info(f"Fetched {len(products)} products for GO term {self.id} using {response_url}")
        Cacher.store_data("go", data_key, products)
        return products

    """
    async def fetch_products_async(self, source):
        #Warning: this function doesn't work if source is GO Annotations File
        if not isinstance(source, GOApi):
            logger.warning("Cannot connect asynchronously, if GOApi is not set as source.")
        async with aiohttp.ClientSession() as session:
            url = source.get_products(self.id)
            response = await session.get(url)
            # TODO: MOVE THIS TO ANNOTATIONPROCESSOR
    """

    def compare_products_to_list(self, products_comparison_list: list):
        """
        Compares self.products (src) to supplied products_comparison_list (ref).
        Returns a list of lists:
          - [0]: a list of products present in self.goterms but not in products_comparison_list
          - [1]: a list of products present in products_comparison_list and not in self.goterms
        """
        products_in_src_and_not_in_ref = []
        products_in_ref_and_not_in_src = []

        for src_product in self.products:
            if src_product not in products_comparison_list:
                products_in_src_and_not_in_ref.append(src_product)

        for ref_product in products_comparison_list:
            if ref_product not in self.products:
                products_in_ref_and_not_in_src.append(ref_product)

        return [products_in_src_and_not_in_ref, products_in_ref_and_not_in_src]

    def compare_products(self, goterm_comparison):
        """
        Compares self.products (src) to supplied goterm_comparison.products (ref).
        Returns a list of lists:
          - [0]: a list of products present in self.goterms but not in goterm_comparison.products
          - [1]: a list of products present in goterm_comparison.products and not in self.goterms
        """
        if isinstance(goterm_comparison, GOTerm):
            return self.compare_products_to_list(
                products_comparison_list=goterm_comparison.products
            )

    def add_process(self, process: Dict):
        if process["process"] not in self.processes:
            self.processes.append(process)

    @classmethod
    def from_dict(cls, d: dict):
        """
        Creates a GOTerm object from a dictionary.

        Args:
            d (dict): A dictionary containing the GO term data.

        Returns:
            A new instance of the GOTerm class.
        """
        # TODO: loop over variables instead of this hardcoded way !!!! using dir(self)

        goterm = GOTerm(id="")
        for attr_name, attr_value in d.items():
            if hasattr(goterm, attr_name):
                # processes must be a list!
                if attr_name == "processes" and isinstance(attr_value, dict):
                    attr_value = [attr_value]
                # check that weight is not passed as string
                if attr_name == "weight" and isinstance(attr_value, str):
                    attr_value = int(attr_value)
                # set the attribute
                setattr(goterm, attr_name, attr_value)
            else:
                logger.warning(f"GO Term class has no attribute name {attr_name}!")
        return goterm
        # old:
        # goterm = cls(
        #    id=d['id'],
        #    processes=d['processes'],
        #    name=d.get('name'),
        #    description=d.get('description'),
        #    weight=d.get('weight', 1.0),
        #    products=d.get('products', []),
        #    http_error_codes=d.get('http_error_codes', {}),
        #    category=d.get('category',None),
        #    parent_term_ids=d.get('parent_term_ids', None),
        #    is_obsolete=d.get('is_obsolete', False)
        #    )
        # return goterm

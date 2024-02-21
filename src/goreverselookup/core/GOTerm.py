import aiohttp, asyncio
from typing import Set, List, Dict, Optional
import json

from ..web_apis.GOApi import GOApi
from ..parse.GOAFParser import GOAnnotationsFile
from ..util.CacheUtil import Cacher
from .ModelSettings import ModelSettings,OrganismInfo
from .ModelStats import ModelStats

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger

class GOTerm:
    def __init__(self, id: str, SOIs: List[Dict] = None, name: Optional[str] = None, description: Optional[str] = None, category: Optional[str] = None, parent_term_ids: Optional[List[str]] = None, is_obsolete:bool = False, weight: float = 1.0, products: List[str] = [], http_error_codes:dict={}):
        """
        A class representing a Gene Ontology term.

        Args:
            id (str): The ID of the GO term.
            SOIs (List[Dict]): [{"SOI" : "angio", "direction" : "+"}]
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            weight (float): The weight of the GO term.
            products (list): Products associated with the term (optional).
            category (str): biological_process, molecular_activity or cellular_component
            parent_term_ids (list[str]): GO ids of the parent terms (parsed from .obo)
            is_obsolete (bool): if the term is labelled as obsolete in the .obo file
        """
        self.id = id
        self.SOIs = SOIs if isinstance(SOIs, list) else [SOIs]
        self.name = name
        self.description = description
        if weight == None: #bugfix for any mistakes
            weight = 1.0
        self.weight = float(weight)
        self.num_products = 0
        self.products = products
        self.products_taxa_dict = {} # a link between a gene id and a belonging taxon
        self.http_error_codes = {} # used for errors happening during server querying; for example, a dict pair 'products': "HTTP Error ..." signifies an http error when querying for GO Term's products
        self.category = category
        self.parent_term_ids = parent_term_ids
        self.is_obsolete = is_obsolete

    def update(self, goterm, overwrite_existing:bool = False):
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
            assert isinstance(goterm,GOTerm)
            for attr_name in dir(goterm): # loop through all class attributes of input goterm
                if not callable(getattr(goterm,attr_name)) and not attr_name.startswith("__"):
                    if hasattr(self, attr_name): # check if the same attribute exists in self
                        if overwrite_existing == True:
                            setattr(self, getattr(goterm,attr_name)) # update self attribute with goterm's attribute
                        else: # overwrite_existing == False
                            attr_value_self = getattr(self, attr_name)
                            if attr_value_self == None or attr_value_self == False:
                                setattr(self, attr_name, getattr(goterm,attr_name)) # update self attribute with goterm's attribute only if self attribute is None or False

        assert isinstance(goterm, GOTerm)
        if self.id == goterm.id: # perform update only if ids match
            overwrite_attributes(goterm=goterm)
        else:
            raise Exception(f"Input GOTerm id {goterm.id} doesn't match with current instance GOTerm id {self.id}")
        
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
        return GOTerm(self.id, self.SOIs, self.name, self.description, self.weight, self.products, self.http_error_codes)
    
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
                self.name = data['label']
            if "definition" in data:
                self.description = data['definition']
            logger.info(f"Fetched name and description for GO term {self.id}")

    async def fetch_name_description_async(self, api: GOApi, session:aiohttp.ClientSession, req_delay=0.1):
        url = api.get_data(self.id, get_url_only=True)
        response = await session.get(url)
        if response.status == 200:
            response_content = await response.read() # response.read() this ensures that response content is fully read before attempting to parse it as JSON
            data = json.loads(response_content)
            # data = await response_content.json()
            await asyncio.sleep(req_delay)
            if "label" in data:
                self.name = data['label']
            if "definition" in data:
                self.description = data['definition']
            # logger.info(f"Fetched name and description for GO term {self.id}")
            # print out only 15 desc chars not to clutter console
            if self.description is not None:
                if len(self.description) > 15:
                    logger.info(f"{self.id}: name = {self.name}, description = {self.description[:15]}...")
            elif self.id is not None and self.name is not None:
                logger.info(f"{self.id}: name = {self.name}, description = None ...")
        else:
            logger.info(f"Query for url {url} failed with response code {response.status}")
    
    def fetch_products(self, source, model_settings:ModelSettings):
        """
        Fetches UniProtKB products associated with a GO Term and sets the "products" member field of the GO Term to a list of all associated products.
        The link used to query for the response is http://api.geneontology.org/api/bioentity/function/{term_id}/genes.

        Parameters:
          - source: can either be a GOApi instance (web-based download) or a GOAnnotationFile isntance (file-based download)
          - (ModelSettings) model_settings: the model settings of this model, used for parsing of target_organism and ortholog organisms
        
        Usage and calling:
            source = GOApi() OR source = GOAnnotationFile()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_products(source)
        """
        if isinstance(source, GOApi):
            r = source.get_products(self.id, model_settings=model_settings)
            if r is not None:
                products = r[0]
                products_taxa_dict = r[1]
                self.products = products
                self.products_taxa_dict = products_taxa_dict
            else:
                self.products = None
                self.products_taxa_dict = None

        elif isinstance(source, GOAnnotationsFile):
            # TODO: reimplement or remove this
            products = source.get_products(self.id)
            if products:
                self.products = products
    
    async def fetch_products_async_v1(self, api:GOApi, delay = 0.0):
        """
        Asynchronously queries the products using api.get_products_async, using the
        await keyword to wait for the request to finish inside the api.
        
        delay is in seconds
        """
        await asyncio.sleep(delay)
        # products = await api.get_products_async(self.id)
        products = await api.get_products_async_notimeout(self.id) # testing variant of the function, comment this later
        if isinstance(products,list) and products != []:
            self.products = products
        elif isinstance(products,str):
            # str reports an HTTP Error code
            if "HTTP Error" in products:
                error_report = products
                self.http_error_codes["products"] = error_report
    
    async def fetch_products_async_v3(self, session:aiohttp.ClientSession, model_settings:ModelSettings, request_params={"rows": 10000000}, req_delay=0.5, max_retries = 3):
        """
        A better variant of get_products_async. Doesn't include timeout in the url request, no retries.
        Doesn't create own ClientSession, but relies on external ClientSession, hence doesn't overload the server as does the get_products_async_notimeout function.
        
        Warning: DO NOT CHANGE request_params={"rows": 10000000}. Decrementing the rows WILL lead to fewer annotations being queried without raising any exceptions in code!

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
        if model_settings.target_organism is None:
            raise Exception("Target organism was not specified in input.txt. Make sure to specify a target organism in the 'settings' section of input.txt!")
        if request_params is not None:
            if 'rows' in request_params:
                if request_params['rows'] < 10000000:
                    logger.warning(f"MAJOR WARNING: Rows specified in request params ({request_params['rows']}) are less than {10000000}. You risk missing out important anootations!")
            if 'rows' not in request_params:
                logger.warning(f"MAJOR WARNING! You did not specify 'rows' in request params. You risk missing out important annotations!")
        
        # data key is in the format [class_name][function_name][function_params]
        data_key = f"[{self.__class__.__name__}][{self.fetch_products_async_v3.__name__}][go_id={self.id}][target_organism={model_settings.target_organism.ncbi_id_full}][orthologs={model_settings.ortholog_organisms_ncbi_full_ids}]"
        
        previous_data_taxa_dict = Cacher.get_data("go", f"{data_key}_products-taxa-dict")
        self.products_taxa_dict = previous_data_taxa_dict
        
        previous_data = Cacher.get_data("go", data_key)
        if previous_data != None:
            logger.debug(f"Cached previous product fetch data for {self.id}")
            ModelStats.goterm_product_query_results[self.id] = previous_data
            self.products = previous_data
            self.products.sort()
            self.num_products = len(self.products)
            return previous_data

        approved_dbs_and_taxa = {} # databases are keys, taxon ids are associated lists
        approved_dbs_and_taxa['UniProtKB'] = [] # create uniprotkb by default
        
        # add target organism
        if model_settings.target_organism.database in approved_dbs_and_taxa:
            # existing database key -> add taxon!
            if model_settings.target_organism.ncbi_id_full not in approved_dbs_and_taxa[model_settings.target_organism.database]:
                approved_dbs_and_taxa[model_settings.target_organism.database] += [model_settings.target_organism.ncbi_id_full]
        else:
            # new database key
            approved_dbs_and_taxa[model_settings.target_organism.database] = [model_settings.target_organism.ncbi_id_full]
        
        # add ortholog organisms to approved_dbs_and_taxa
        if model_settings.ortholog_organisms is not None:
            for ortholog_organism_id, ortholog_organism_object in model_settings.ortholog_organisms.items():
                assert isinstance(ortholog_organism_object, OrganismInfo)
                ortholog_organism = ortholog_organism_object
                if ortholog_organism.database in approved_dbs_and_taxa:
                    # existing database key -> add taxon
                    if ortholog_organism.ncbi_id_full not in approved_dbs_and_taxa[ortholog_organism.database]:
                        approved_dbs_and_taxa[ortholog_organism.database] += [ortholog_organism.ncbi_id_full]
                
                    # automatically add ortholog organism taxon to the UniProtKB database query
                    if ortholog_organism.ncbi_id_full not in approved_dbs_and_taxa['UniProtKB']:
                        approved_dbs_and_taxa['UniProtKB'] += [ortholog_organism.ncbi_id_full]
                else:
                    approved_dbs_and_taxa[ortholog_organism.database] = [ortholog_organism.ncbi_id_full]
        
        url = f"http://api.geneontology.org/api/bioentity/function/{self.id}/genes"
        params = request_params # 10k rows resulted in 56 mismatches for querying products for 200 goterms (compared to reference model, loaded from synchronous query data)
        
        retries = 0
        data = None
        for _ in range(model_settings.goterm_gene_query_max_retries):
            possible_http_error_text = ""
            if retries == (max_retries-1):
                logger.warning(f"Exceeded max retries when parsing {self.id}")
                ModelStats.goterm_product_query_results[self.id] = f"Error: Exceeded max retries. Error info: {possible_http_error_text}"
                return None
            retries +=1

            previous_response = Cacher.get_data("url", url)
            if previous_response != None:
                data = previous_response
                break # previous response json was cached, break the loop
            else:
                try:
                    await asyncio.sleep(req_delay)
                    response = await session.get(url, params=params, timeout=model_settings.goterm_gene_query_timeout)
                    if response.status != 200: # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                        possible_http_error_text = f"HTTP Error when parsing {self.id}. Response status = {response.status}"
                        logger.warning(possible_http_error_text)
                        continue
                    # data = await response.json()
                    response_content = await response.read()
                    data = json.loads(response_content)
                    if data != None:
                        Cacher.store_data("url", url, data)
                        logger.debug(f"Cached async product fetch data for {self.id}")
                        break # reponse json was obtained, break the loop
                except (aiohttp.ClientConnectionError, aiohttp.ClientPayloadError, asyncio.TimeoutError) as e:
                    logger.warning(f"Error when fetching products for {self.id}: {type(e).__name__}")
                    logger.warning(f"  - attempted url: {url}")
                    possible_http_error_text = f"{e}"
            
        products_set = set()
        products_taxa_dict = {}
        _d_unique_dbs = set() # unique databases of associations; eg. list of all unique assoc['subject']['id']
        
        associations = data.get('associations', []) if data is not None else []
        
        for assoc in associations:
            # if evidence is not confirmed, continue to next iteration
            evidence_confirmed, evidence_code_eco_id = GOApi.check_GO_association_evidence_code_validity(assoc, model_settings.valid_evidence_codes)
            if not evidence_confirmed:
                continue
   
            _d_unique_dbs.add(assoc['subject']['id'].split(":")[0])
            _d_db = assoc['subject']['id'].split(":")[0]
            _d_taxon = assoc['subject']['taxon']['id']

            if assoc['object']['id'] == self.id:
                for database,taxa in approved_dbs_and_taxa.items():
                    if database in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in taxa):
                        product_id = assoc['subject']['id']
                        products_set.add(product_id)
                        products_taxa_dict[product_id] = assoc['subject']['taxon']['id']
        
        products = list(products_set)
        if products == []:
            logger.warning(f"Found no products for GO Term {self.id} (name = {self.name})!")
            #if len(data) < 500:
            #    logger.debug(f"Response json: {data}")

        self.products = products
        self.products.sort() # alphabetically order
        self.num_products = len(products)
        self.products_taxa_dict = products_taxa_dict
        # logger.debug(f"Active session connections: {len(session.connector._conns)}")
        logger.info(f"Fetched {len(self.products)} products for GO term {self.id} from {len(_d_unique_dbs)} unique databases ({_d_unique_dbs})")
        Cacher.store_data("go", data_key, products)
        Cacher.store_data("go", f"{data_key}_products-taxa-dict", self.products_taxa_dict)
        ModelStats.goterm_product_query_results[self.id] = self.products
        return products

    def compare_products_to_list(self, products_comparison_list:list):
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
        
        return[products_in_src_and_not_in_ref, products_in_ref_and_not_in_src]


    def compare_products(self, goterm_comparison):
        """
        Compares self.products (src) to supplied goterm_comparison.products (ref).
        Returns a list of lists:
          - [0]: a list of products present in self.goterms but not in goterm_comparison.products
          - [1]: a list of products present in goterm_comparison.products and not in self.goterms
        """
        if isinstance(goterm_comparison, GOTerm):
            return self.compare_products_to_list(products_comparison_list = goterm_comparison.products)

    def add_SOI(self, SOI: Dict):
        if not SOI['SOI'] in self.SOIs:
            self.SOIs.append(SOI)

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
        for attr_name,attr_value in d.items():
            if hasattr(goterm, attr_name):
                if attr_name == 'SOIs' and isinstance(attr_value,dict):
                    attr_value = [attr_value]
                # check that weight is not passed as string
                if attr_name == 'weight' and isinstance(attr_value,str):
                    attr_value = int(attr_value)
                # set the attribute
                setattr(goterm, attr_name, attr_value)
            else:
                logger.warning(f"GO Term class has no attribute name {attr_name}!")
        return goterm
    
    def to_json(self):
        json_data = {}
        for attr_name, attr_value in vars(self).items():
            if attr_value is None:
                # this happens for example when only attr_name is defined without a value; e.g. "ortholog_organisms" without any following ortholog organisms.
                continue
            # custom handling for target_organism and ortholog_organisms, as they are code objects -> convert them to json
            if not callable(attr_value) and not attr_name.startswith("__"):
                # append to json_data result dict
                json_data[attr_name] = attr_value
        return json_data
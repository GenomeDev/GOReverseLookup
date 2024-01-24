import requests
from requests.adapters import HTTPAdapter, Retry
from json import JSONDecodeError
import time
import aiohttp
import asyncio
import json

from ..util.FileUtil import FileUtil
from ..core.ModelStats import ModelStats
from ..util.CacheUtil import Cacher
from ..core.ModelSettings import ModelSettings,OrganismInfo

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class GOApi:
    """
    This class enables the user to interact with the Gene Ontology database via http requests.
    """

    def __init__(self):
        # Set up a retrying session
        retry_strategy = Retry(
            total=3, 
            status_forcelist=[429, 500, 502, 503, 504], 
            backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session
        self.api_version = None
    
    def get_GO_version(self):
        """
        Returns the GO version. This sets the self.api_version parameter.
        """
        url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/about"
        r = self.s.get(url, headers={ "Accept" : "application/json"})
        response_json = r.json()
        return response_json['go']['timestamp']
            

    def get_data(self, term_id, get_url_only=False):
        """
        Fetches term data for a given term ID from the Gene Ontology API using http://api.geneontology.org/api/ontology/term/{term_id},
        example of a term_id is GO:1903589.

        If get_url_only == True, this will only return the url.

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
        url = f"http://api.geneontology.org/api/ontology/term/{term_id}"
        params = {}
        if get_url_only:
            return url
        logger.debug(f"Querying: {url}")
        try:
            response = self.s.get(url, params=params, timeout=5)
            if response.ok:
                data = response.json()
                return data
            else:
                logger.warning(f"Error: {response.status_code} - {response.reason}")
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error: {e}")
            return None

    def get_products(
        self, 
        term_id, 
        model_settings:ModelSettings, 
        get_url_only=False, 
        get_response_only=False, 
        request_params={"rows": 10000000}
    ):
        """
        Fetches product IDs (gene ids) associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        The request uses this link: http://api.geneontology.org/api/bioentity/function/{term_id}/genes

        Returns: a list of two elements
          - [0]: (string as json) data: a json string, representing the api request response
          - [1]: (dict) products_taxa_dict
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
        data_key = f"[{self.__class__.__name__}][{self.get_products.__name__}][go_id={term_id}][target_organism={model_settings.target_organism.ncbi_id_full}][orthologs={model_settings.ortholog_organisms_ncbi_full_ids}]"
        prev_data = Cacher.get_data("go", data_key=data_key)
        previous_data_taxa_dict = Cacher.get_data("go", f"{data_key}_products-taxa-dict")

        if prev_data is not None and previous_data_taxa_dict is not None:
            logger.debug(f"Found cached previous product fetch data for {term_id}")
            ModelStats.goterm_product_query_results[term_id] = prev_data
            return [prev_data, previous_data_taxa_dict]

        approved_dbs_and_taxa = {} # databases are keys, taxon ids are associated lists
        approved_dbs_and_taxa['UniProtKB'] = [] # create uniprotkb by default
        products_taxa_dict = {}

        # add target organism
        if model_settings.target_organism.database in approved_dbs_and_taxa:
            # existing database key -> add taxon!
            if model_settings.target_organism.ncbi_id_full not in approved_dbs_and_taxa[model_settings.target_organism.database]:
                approved_dbs_and_taxa[model_settings.target_organism.database] += [model_settings.target_organism.ncbi_id_full]
        else:
            # new database key
            approved_dbs_and_taxa[model_settings.target_organism.database] = [model_settings.target_organism.ncbi_id_full]
        
        # add ortholog organisms
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
         
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = request_params

        # used in async
        if get_url_only is True:
            # create a request object with the base url and params
            request = requests.Request("GET", url, params=params)
            # prepare the request
            prepared_request = self.s.prepare_request(request)
            # get the fully constructed url with parameters
            url = prepared_request.url
            return url

        max_retries = 5  # try api requests for max 5 times
        for i in range(model_settings.goterm_gene_query_max_retries):
            try:
                response = self.s.get(url, params=params, timeout=model_settings.goterm_gene_query_timeout)
                response.raise_for_status()

                data = response.json()
                if get_response_only == True:
                    return data
                
                products_set = set()
                _d_unique_dbs = set() # unique databases of associations; eg. list of all unique assoc['subject']['id']
        
                for assoc in data['associations']:
                    evidence_confirmed, evidence_code_eco_id = GOApi.check_GO_association_evidence_code_validity(assoc, model_settings.valid_evidence_codes)
                    if not evidence_confirmed:
                        continue

                    _d_unique_dbs.add(assoc['subject']['id'].split(":")[0])
                    _d_db = assoc['subject']['id'].split(":")[0]
                    _d_taxon = assoc['subject']['taxon']['id']

                    if assoc['object']['id'] in term_id or term_id in assoc['object']['id']:
                        for database,taxa in approved_dbs_and_taxa.items():
                            if database in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in taxa):
                                product_id = assoc['subject']['id']
                                products_set.add(product_id)
                                products_taxa_dict[product_id] = assoc['subject']['taxon']['id']
                
                products = list(products_set)
                logger.info(f"{term_id}: fetched {len(products)} products.")
                Cacher.store_data(data_location="go", data_key=data_key, data_value=products)
                Cacher.store_data("go", f"{data_key}_products-taxa-dict", products_taxa_dict)
                ModelStats.goterm_product_query_results[term_id] = products
                return [products, products_taxa_dict]
            
            except (requests.exceptions.RequestException, JSONDecodeError) as e:
                if i == (max_retries - 1):  # this was the last http request, it failed
                    logger.error(f"Experienced an http exception or a JSONDecodeError while fetching products for {term_id}")
                    ModelStats.goterm_product_query_results[term_id] = f"Error: {type(e).__name__}"
                    
                    error_log_filepath = FileUtil.find_win_abs_filepath("log_output/error_log")
                    error_type = type(e).__name__
                    error_text = str(e)

                    logger.error(f"Exception type: {error_type}")
                    logger.error(f"Exception text: {error_text}")
                    logger.error(f"Debug report was written to: {error_log_filepath}")

                    with open(error_log_filepath, "a+") as f:
                        f.write(f"Fetch products error for: {term_id}\n")
                        f.write(f"Exception: {error_type}\n")
                        f.write(f"Cause: {error_text}\n")
                        f.write("\n\n\n")
                        f.write("------------------------------\n")
                else:
                    # time.sleep(500) # sleep 500ms before trying another http request
                    time.sleep(0.5)  # time.sleep is in SECONDS !!!
                return None
            
        logger.warning(f"Exceeded max retries when querying products for {term_id}")
        ModelStats.goterm_product_query_results[term_id] = f"Error: Exceeded max retries."
        return None

    async def get_products_async(self, term_id, model_settings:ModelSettings):
        """
        Fetches product IDs associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        This function works asynchronously, much faster than it's synchronous 'get_products' counterpart.

        The request uses this link: http://api.geneontology.org/api/bioentity/function/{term_id}/genes

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
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
        
        # add ortholog organisms
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
        
        MAX_RETRIES = 5
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = {"rows": 100000}

        global request_iterations
        request_iterations = 0  # global variable request_iterations to keep track of the amount of requests submitted to the server (maximum is MAX_RETRIES); a harsh bugfix

        # as per: https://stackoverflow.com/questions/51248714/aiohttp-client-exception-serverdisconnectederror-is-this-the-api-servers-issu
        connector = aiohttp.TCPConnector(limit=20)  # default limit is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            # for i in range(MAX_RETRIES):
            # while i < MAX_RETRIES: # due to the async nature, each iteration resets i; hence "i" is useless -> bugfix: global variable request_iterations
            while request_iterations < MAX_RETRIES:
                try:
                    request_iterations += 1
                    response = await session.get(url, params=params, timeout=7)
                    response.raise_for_status()  # checks for anything other than status 200
                    # data = await response.json()
                    response_content = await response.read()
                    data = json.loads(response_content)

                    products_set = set()
                    _d_unique_dbs = set() # unique databases of associations; eg. list of all unique assoc['subject']['id']
        
                    for assoc in data['associations']:
                        _d_unique_dbs.add(assoc['subject']['id'].split(":")[0])
                        _d_db = assoc['subject']['id'].split(":")[0]
                        _d_taxon = assoc['subject']['taxon']['id']
                        if "UniProtKB" in _d_db and "9606" not in _d_taxon: # for debug purposes
                            # logger.info(f"UniProtKB - {_d_taxon}")
                            pass
                        if assoc['object']['id'] == self.id:
                            for database,taxa in approved_dbs_and_taxa.items():
                                if database in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in taxa):
                                    product_id = assoc['subject']['id']
                                    products_set.add(product_id)

                    products = list(products_set)
                    logger.info(f"Fetched products for GO term {term_id}")
                    request_iterations = 0  # reset
                    return products
                except (
                    requests.exceptions.RequestException,
                    JSONDecodeError,
                    asyncio.exceptions.TimeoutError,
                    aiohttp.ClientResponseError,
                ) as e:
                    # logger.error(f"TimoutError on retry attempt {request_iterations}. Exception: {e}")
                    # i += 1
                    # if i == (MAX_RETRIES - 1): # this was the last http request, it failed
                    # if request_iterations == (MAX_RETRIES - 1):
                    if (
                        request_iterations == MAX_RETRIES
                    ):  # due to while loop logic we don't substract 1
                        error_log_filepath = FileUtil.find_win_abs_filepath(
                            "log_output/error_log"
                        )
                        error_type = type(e).__name__
                        error_text = str(e)

                        # logger.error(f"Exception type: {error_type}")
                        # logger.error(f"Exception text: {error_text}")
                        # logger.error(f"Debug report was written to: {error_log_filepath}")
                        logger.error(
                            f"https error for {term_id}, error_type = {error_type},"
                            f" error_text = {error_text}"
                        )

                        with open(error_log_filepath, "a+") as f:
                            f.write(f"Fetch products error for: {term_id}\n")
                            f.write(f"Exception: {error_type}\n")
                            f.write(f"Cause: {error_text}\n")
                            f.write("\n\n\n")
                            f.write("------------------------------\n")
                    else:
                        # time.sleep(0.5)
                        time.sleep(1)  # maybe with 1s the server won't start to block?
            # reset
            request_iterations = 0

    async def get_products_async_notimeout(self, term_id):
        """
        A testing variant of get_products_async. Doesn't include timeout in the url request, no retries.
        """
        APPROVED_DATABASES = [
            ["UniProtKB", ["NCBITaxon:9606"]],
            ["ZFIN", ["NCBITaxon:7955"]],
            # ["RNAcentral", ["NCBITaxon:9606"]],
            ["Xenbase", ["NCBITaxon:8364"]],
            ["MGI", ["NCBITaxon:10090"]],
            ["RGD", ["NCBITaxon:10116"]],
        ]
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = {
            "rows": 20000
        }  # 10k rows resulted in 56 mismatches for querying products for 200 goterms (compared to reference model, loaded from synchronous query data)
        # DELAY = 1 # 1 second delay between requests

        # as per: https://stackoverflow.com/questions/51248714/aiohttp-client-exception-serverdisconnectederror-is-this-the-api-servers-issu
        connector = aiohttp.TCPConnector(
            limit=20, limit_per_host=20
        )  # default limit is 100
        # as per: https://stackoverflow.com/questions/64534844/python-asyncio-aiohttp-timeout; DOESNT WORK!
        # session_timeout =   aiohttp.ClientTimeout(total=None,sock_connect=10,sock_read=10) -> async with aiohttp.ClientSession(connector=connector, timeout=session_timeout) as session;
        # https://github.com/aio-libs/aiohttp/issues/3187 -> 504 gateways are server-limited !

        ### POSSIBLE ERROR SOLUTION ### [TODO: continue from here]
        # Current algorithm creates one aiohttp.ClientSession FOR EACH GOTERM. Therefore, each ClientSession only has one connection,
        # and the checks for connection limiting aren't enforeced. During runtime, there can be as many as 200 (as many as there are goterms)
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

        async with aiohttp.ClientSession(connector=connector) as session:
            response = await session.get(url, params=params)
            # response.raise_for_status() # checks for anything other than status 200
            if response.status != 200:  # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                logger.warning(f"HTTP Error when parsing {term_id}. Response status ={response.status}")
                return (f"HTTP Error: status = {response.status}, reason = {response.reason}")

            # data = await response.json()
            response_content = await response.read()
            data = json.loads(response_content)
            products_set = set()
            for assoc in data["associations"]:
                if assoc["object"]["id"] == term_id and any(
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
            logger.info(f"Fetched products for GO term {term_id}")
            return products

    def get_goterms(
        self,
        gene_id: str,
        go_categories: list = [
            "molecular_activity",
            "biological_process",
            "cellular_component",
        ],
        approved_taxa=["NCBITaxon:9696"],
        request_params={"rows": 10000000},
        model_settings:ModelSettings=None
    ):
        """
        Gets all GO Terms associated with 'gene_id' in the form of a list.

        Parameters:
          - (str) gene_id: The full gene id (eg. UniProtKB:P15692)
          - (list) go_categories: a list of valid categories. All possible categories are 'molecular_activity', 'biological_process', 'cellular_component'.
                                  All categories are accepted by default.
          - () request_params: leave it be. Shortening may cause incomplete JSON objects to be returned.
          - (list) approved_taxa: All the taxa that can be returned. If get_goterms is used inside fisher_exact_test for Fisher scoring (ModelSettings.fisher_test_use_online_query == True),
                                  then the value of this parameter greatly determines the amount of GO Terms (associated to a gene) that are returned. Specifically, this parameter determines the
                                  num_goterms_product_general value of the Fisher exact test contingency table. Only include the taxon (or taxa) which is (are) part of the research. If you are interested
                                  in statistically significant genes for Homo Sapiens, then only include the Homo Sapiens NCBI Taxon.
                                  The taxon (taxa) should be in the form of a list, each taxon should be a full NCBITaxon, such as: ["NCBITaxon:9696"]
          - (ModelSettings) model_settings: needed for valid reference codes
          
        To carry out the query request, the following url is used:
            http://api.geneontology.org/api/bioentity/gene/{gene_id}/function
        """
        url = f"http://api.geneontology.org/api/bioentity/gene/{gene_id}/function"
        response = requests.get(url, params=request_params)
        result_go_terms = []

        if response.status_code == 200:
            response_json = response.json()
            for assoc in response_json["associations"]:
                # if evidence is not confirmed, continue to next iteration
                evidence_confirmed, evidence_code_eco_id = GOApi.check_GO_association_evidence_code_validity(assoc, model_settings.valid_evidence_codes)
                if not evidence_confirmed:
                    continue
                if assoc["subject"]["taxon"]["id"] in approved_taxa:
                    if assoc["object"]["category"][0] in go_categories:
                        go_id = assoc["object"]["id"]
                        if go_id is not None:
                            result_go_terms.append(go_id)
            return result_go_terms
        else:
            logger.warning(f"Response error when querying GO Terms for {gene_id}!")
            return None
    
    @classmethod
    def check_GO_association_evidence_code_validity(cls, assoc, valid_evidence_codes):
        """
        Checks if an association (between a GO term and a gene) has a matching evidence code as are the
        specified valid_evidence_codes list.

        Parameters:
          - assoc: dictionary representing the association. Associations are obtained for example using http://api.geneontology.org/api/bioentity/function/{term_id}/genes, querying for response.json()
                    and then using response_json["associations"]
          - valid_evidence_codes: a list of valid evidence codes. They must be in the ECO format. It is preferred to instantiate ModelSettings with evidence codes and pass evidence codes from ModelSettings.
        
        Returns: list
          - [0]: True or False, if the association has a valid evidence code
          - [1]: the evidence code of the association
        """
        evidence_code_eco_id = ""
        evidence_confirmed = False # if any of the evidence codes for this object are among (ReverseLookup).model_settings.valid_evidence_codes
        evidence_code = assoc.get('evidence', None)
        if evidence_code is not None:
            if isinstance(evidence_code, list):
                for ev_c in evidence_code:
                    evidence_code_eco_id = ev_c
                    if ev_c in valid_evidence_codes:
                        evidence_confirmed = True
                        break
            elif isinstance(evidence_code, str):
                evidence_code_eco_id = evidence_code
                if evidence_code in valid_evidence_codes:
                    evidence_confirmed = True
        else: # evidence_code is None -> check evidence types
            evidence_types = assoc.get('evidence_types')
            if evidence_types is not None:
                for evidence_type in evidence_types:
                    evidence_code_eco_id = evidence_type.get('id', None) # get ECO id
                    if evidence_code_eco_id in valid_evidence_codes:
                        evidence_confirmed = True
                        break
        return [evidence_confirmed, evidence_code_eco_id]

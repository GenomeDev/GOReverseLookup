import requests
from requests.adapters import HTTPAdapter, Retry
from json import JSONDecodeError
import time
import aiohttp
import asyncio

from ..util.FileUtil import FileUtil

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class GOApi:
    """
    This class enables the user to interact with the Gene Ontology database via http requests.
    """

    def __init__(self):
        # Set up a retrying session
        retry_strategy = Retry(
            total=3, status_forcelist=[429, 500, 502, 503, 504], backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session

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
        self, term_id, get_url_only=False, get_response_only=False, request_params={"rows": 10000000}
    ):
        """
        Fetches product IDs (gene ids) associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        The request uses this link: http://api.geneontology.org/api/bioentity/function/{term_id}/genes

        Returns:
          - (string as json) data: a json string, representing the api request response
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

        products_set = set()
        max_retries = 5  # try api requests for max 5 times
        for i in range(max_retries):
            try:
                response = self.s.get(url, params=params, timeout=5)
                response.raise_for_status()

                json = response.json()
                if get_response_only == True:
                    return json
                
                for assoc in json["associations"]:
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

            except (requests.exceptions.RequestException, JSONDecodeError) as e:
                if i == (max_retries - 1):  # this was the last http request, it failed
                    logger.error(
                        "Experienced an http exception or a JSONDecodeError while"
                        f" fetching products for {term_id}"
                    )
                    error_log_filepath = FileUtil.find_win_abs_filepath(
                        "log_output/error_log"
                    )
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

    async def get_products_async(self, term_id):
        """
        Fetches product IDs associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        This function works asynchronously, much faster than it's synchronous 'get_products' counterpart.

        The request uses this link: http://api.geneontology.org/api/bioentity/function/{term_id}/genes

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
        APPROVED_DATABASES = [
            ["UniProtKB", ["NCBITaxon:9606"]],
            ["ZFIN", ["NCBITaxon:7955"]],
            # ["RNAcentral", ["NCBITaxon:9606"]],
            ["Xenbase", ["NCBITaxon:8364"]],
            ["MGI", ["NCBITaxon:10090"]],
            ["RGD", ["NCBITaxon:10116"]],
        ]
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
                    data = await response.json()
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
            if (
                response.status != 200
            ):  # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                logger.warning(
                    f"HTTP Error when parsing {term_id}. Response status ="
                    f" {response.status}"
                )
                return (
                    f"HTTP Error: status = {response.status}, reason ="
                    f" {response.reason}"
                )

            data = await response.json()
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

        To carry out the query request, the following url is used:
            http://api.geneontology.org/api/bioentity/gene/{gene_id}/function
        """
        url = f"http://api.geneontology.org/api/bioentity/gene/{gene_id}/function"
        response = requests.get(url, params=request_params)
        result_go_terms = []

        if response.status_code == 200:
            response_json = response.json()
            for assoc in response_json["associations"]:
                if assoc["subject"]["taxon"]["id"] in approved_taxa:
                    if assoc["object"]["category"][0] in go_categories:
                        go_id = assoc["object"]["id"]
                        if go_id is not None:
                            result_go_terms.append(go_id)
            return result_go_terms
        else:
            logger.warning(f"Response error when querying GO Terms for {gene_id}!")
            return None

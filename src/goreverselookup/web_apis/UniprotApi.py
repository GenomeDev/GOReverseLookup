import requests
from requests.adapters import HTTPAdapter, Retry
import aiohttp
import asyncio

from ..util.CacheUtil import Cacher

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class UniProtApi:
    """
    This class enables the user to interact with the UniProtKB database via http requests.
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
        self.uniprot_query_exceptions = []
        self.async_request_sleep_delay = 0.5

    def get_uniprot_id(self, gene_name, get_url_only=False):
        """
        Given a gene name, returns the corresponding UniProt ID using the UniProt API.

        Parameters:
        - gene_name (str): name of the gene to search for.
        - TODO: retries (int): maximum number of times to retry the request in case of network errors.
        - TODO: timeout (int): timeout in seconds for the request.
        - get_url_only: only return the query url without performing the query

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections

        Returns:
        - str: UniProt ID if found, None otherwise OR the query url, if get_url_only is True
        """
        # data key is in the format [class_name][function_name][function_params]
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_id.__name__}][gene_name={gene_name}]"
        previous_uniprot_id = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_uniprot_id is not None:
            logger.debug(
                f"Cached uniprot id {previous_uniprot_id} for gene name {gene_name}"
            )
            return previous_uniprot_id

        # Define the URL to query the UniProt API
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl"

        if get_url_only:
            return url

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            # Try the request up to `retries` times
            try:
                # Make the request and raise an exception if the response status is not 200 OK
                response = self.s.get(url, timeout=5)
                response.raise_for_status()
                response_json = response.json()
                Cacher.store_data("url", url, response_json)
            except requests.exceptions.RequestException:
                # if there was an error with the HTTP request, log a warning
                logger.warning(f"Failed to fetch UniProt data for {gene_name}")
                return None

        # Parse the response JSON and get the list of results
        # results = response.json()["results"]
        results = response_json["results"]

        # If no results were found, return None
        if len(results) == 0:
            return None

        # If only one result was found, accept it automatically
        elif len(results) == 1:
            uniprot_id = results[0]["primaryAccession"]
            logger.info(
                f"Auto accepted {gene_name} -> {uniprot_id}. Reason: Only 1 result."
            )
            return_value = "UniProtKB:" + uniprot_id
            Cacher.store_data("uniprot", uniprot_data_key, return_value)
            return return_value

        # If multiple results were found, filter out the non-reviewed ones
        reviewed_ids = []
        for result in results:
            # Skip the result if the gene name is not a match
            if gene_name not in result["genes"][0]["geneName"]["value"]:
                continue
            # Skip the result if it is not reviewed
            if "TrEMBL" not in result["entryType"]:
                reviewed_ids.append(result)

        # If no reviewed result was found, return None
        if len(reviewed_ids) == 0:
            return None

        # If only one reviewed result was found, accept it automatically
        elif len(reviewed_ids) == 1:
            uniprot_id = reviewed_ids[0]["primaryAccession"]
            logger.info(
                f"Auto accepted {gene_name} -> {uniprot_id}. Reason: Only 1 reviewed"
                " result."
            )
            return_value = "UniProtKB:" + uniprot_id
            Cacher.store_data("uniprot", uniprot_data_key, return_value)
            return return_value

        # If multiple reviewed results were found, ask the user to choose one
        logger.info(
            f"Multiple reviewed results found for {gene_name}. Please choose the"
            " correct UniProt ID from the following list:"
        )
        for i, result in enumerate(reviewed_ids):
            genes = result["genes"]
            impact_genes = set()
            for gene in genes:
                impact_genes.add(gene["geneName"]["value"])
                if "synonyms" in gene:
                    for synonym in gene["synonyms"]:
                        impact_genes.add(synonym["value"])
            print(f"{i + 1}. {result['primaryAccession']} ({', '.join(impact_genes)})")
        # Get the user's choice and return the corresponding UniProt ID
        # choice = input("> ")  # prompt the user for input, but commented out for now
        choice = "1"  # for testing purposes, use "1" as the user's choice
        if choice.isdigit() and 1 <= int(choice) <= len(
            reviewed_ids
        ):  # check if the user's choice is valid
            # get the UniProt ID of the chosen result and return it
            uniprot_id = reviewed_ids[int(choice) - 1]["primaryAccession"]
            logger.warning(f"Auto-selected first reviewed result for {gene_name}!")
            return_value = "UniProtKB:" + uniprot_id
            Cacher.store_data("uniprot", uniprot_data_key, return_value)
            return return_value
        else:
            # raise an error if the user's choice is not valid
            raise ValueError(f"Invalid choice: {choice}")

    def get_uniprot_info(self, uniprot_id: str) -> dict:
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.

        Parameters:
          - (str) uniprot_id

        If the query is successful, returns the following dictionary:
            {
                "genename": name,
                "description": description,
                "ensg_id": ensg_id,
                "enst_id": enst_id,
                "refseq_nt_id": refseq_nt_id
            }

        This function automatically uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections

        Algorithm:
        This function constructs a uniprot url (https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name)
        and queries the response. A query response example: https://pastebin.com/qrWck3QG. The structure of the query response is:
        {
            "results":[
                {
                    "entryType": "UniProtKB reviewed (Swiss-Prot)",
                    "primaryAccession": "Q9NY91",
                    "organism": {
                        "scientificName": "Homo Sapiens",
                        "commonName": "Human",
                        "taxonId": 9606
                        "lineage": [...]
                    },
                    "proteinDescription": {
                        "recommendedName": {
                            "fullName": {
                                "value": "Probable glucose sensor protein SLC5A4"
                                "evidences": [...]
                            }
                        },
                        "alternativeNames": {
                            "fullName": {
                                ...
                            }
                        }
                    },
                    "genes": [
                        {
                            "genename": {
                                "value": "SLC5A4",
                                "evidences": [...]
                            },
                            "synonyms": [
                                {
                                    "value": "SAAT1"
                                },
                                {
                                    "value": "SGLT3"
                                    "evidences": [...]
                                }
                            ]
                        }
                    ],
                    "uniProtKBCrossReferences": [
                        {
                            "database": "RefSeq",
                            "id": "NP_055042.1",
                            "properties": [
                                {
                                    "key": "NucleotideSequenceId",
                                    "value": "NM_014227.2"
                                },
                                {
                                    "database": "Ensembl",
                                    "id": "ENST00000266086.6",
                                    "properties": [
                                        {
                                            "key": "ProteinId",
                                            "value": "ENSP00000266086.3"
                                        },
                                        {
                                            "key": "GeneId",
                                            "value": "ENSG00000100191.6"
                                        }
                                    ]
                                },
                                {
                                    "database": "MANE-Select",
                                    "id": "ENST00000266086.6",
                                    "properties": [
                                        {
                                            "key": "ProteinId",
                                            "value": "ENSP00000266086.3"
                                        },
                                        {
                                            "key": "RefSeqNucleotideId",
                                            "value": "NM_014227.3"
                                        },
                                        {
                                            "key": "RefSeqProteinId",
                                            "value": "NP_055042.1"
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                },
                {
                # repeat above structure for more results
                }
            ]
        }

        This function then passes the response results (response_json["results"]) to
        _process_uniprot_info_query_results(response_results, uniprot_id), which returns
        the final dictionary with the processed values.
        """
        # Extract UniProt ID if given in "database:identifier" format
        if ":" in uniprot_id:
            uniprot_id = uniprot_id.split(":")[1]

        # Attempt to return previously cached function return value
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_info.__name__}][uniprot_id={uniprot_id}]"
        previous_info = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_info is not None:
            logger.debug(
                f"Returning cached info for uniprot id {uniprot_id}: {previous_info}"
            )
            return previous_info

        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            try:
                response = requests.get(url, timeout=5)
                response.raise_for_status()
                response_json = response.json()
                Cacher.store_data("url", url, response_json)
            except requests.exceptions.RequestException:
                logger.warning(f"Failed to fetch UniProt data for {uniprot_id}")
                return {}

        # TODO: delete this debug!
        logger.debug(f"url = {url}")
        logger.debug(f"response_json = {response_json}")
        results = response_json["results"]
        return_value = self._process_uniprot_info_query_results(results, uniprot_id)

        # cache function result
        if return_value is not None:
            Cacher.store_data("uniprot", uniprot_data_key, return_value)
        return return_value

    def _return_mane_select_values_from_uniprot_query(self, result: dict) -> tuple:
        """
        Given the UniProt search result dictionary, return Ensembl gene ID, Ensembl transcript ID, and RefSeq nucleotide ID for the MANE-select transcript.

        Usage:
        (1) get the uniprot id in question UniProtKB:Q9NY91 -> uniprot_id = Q9NY91
        (2) query info using https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name
        (3) query_result = response_json["results"][0]
        (4) enst_id, refseq_nt_id, ensg_id = _return_mane_select_values_from_uniprot_query(query_result)

        This function is used in the (UniprotAPI).get_info -> _process_uniprot_info_query_results(...) function.
        """
        # uniProtKBCrossReferences structure:
        # "uniProtKBCrossReferences": [
        #        {
        #           "database": "RefSeq",
        #           "id": "NP_055042.1",
        #           "properties": [
        #                {
        #                    "key": "NucleotideSequenceId",
        #                    "value": "NM_014227.2"
        #                }
        #            ]
        #        },
        #        {
        #            "database": "Ensembl",
        #            "id": "ENST00000266086.6",
        #            "properties": [
        #                {
        #                    "key": "ProteinId",
        #                    "value": "ENSP00000266086.3"
        #                },
        #                {
        #                    "key": "GeneId",
        #                    "value": "ENSG00000100191.6"
        #                }
        #            ]
        #        },
        #        {
        #            "database": "MANE-Select",
        #            "id": "ENST00000266086.6",
        #            "properties": [
        #                {
        #                    "key": "ProteinId",
        #                    "value": "ENSP00000266086.3"
        #                },
        #                {
        #                    "key": "RefSeqNucleotideId",
        #                    "value": "NM_014227.3"
        #                },
        #                {
        #                    "key": "RefSeqProteinId",
        #                    "value": "NP_055042.1"
        #                }
        #            ]
        #        }
        #    ]
        # }

        # inside uniProtKBCrossReferences dictionary, find the index of the MANE-Select element. In the above example, the MANE-Select element is the third in array, hence it has index 2 -> [2]
        mane_indices = [
            index
            for (index, d) in enumerate(result["uniProtKBCrossReferences"])
            if d["database"] == "MANE-Select"
        ]
        if (
            len(mane_indices) == 1
        ):  # only one MANE-Select element in uniProtKBCrossReferences
            i = mane_indices[0]
            enst_id = result["uniProtKBCrossReferences"][i]["id"]
            refseq_nt_id = next(
                (
                    entry["value"]
                    for entry in result["uniProtKBCrossReferences"][i]["properties"]
                    if entry["key"] == "RefSeqNucleotideId"
                ),
                None,
            )
            ensg_id = next(
                (
                    next(
                        (
                            sub["value"]
                            for sub in entry["properties"]
                            if sub["key"] == "GeneId"
                        ),
                        None,
                    )
                    for entry in result["uniProtKBCrossReferences"]
                    if (entry["database"] == "Ensembl" and entry["id"] == enst_id)
                ),
                None,
            )
            return ensg_id, enst_id, refseq_nt_id
        else:
            return None, None, None

    def _process_uniprot_info_query_results(
        self, results: str, uniprot_id: str
    ) -> dict:
        """
        Processes the results obtained from the get_uniprot_info query.

        Returns the following dictionary:
            {
                "genename": name,
                "description": description,
                "ensg_id": ensg_id,
                "enst_id": enst_id,
                "refseq_nt_id": refseq_nt_id
            }

        See the JSON structure in the comment of the get_uniprot_info function.
        """
        if len(results) == 0:
            return {}
        else:
            # Get values from the UniProt search result
            result = next(
                (entry for entry in results if entry["primaryAccession"] == uniprot_id),
                None,
            )
            name = result["genes"][0]["geneName"]["value"]  # gene name
            if (
                "proteinDescription" in result
                and "recommendedName" in result["proteinDescription"]
                and "fullName" in result["proteinDescription"]["recommendedName"]
                and "value"
                in result["proteinDescription"]["recommendedName"]["fullName"]
            ):
                description = result["proteinDescription"]["recommendedName"][
                    "fullName"
                ]["value"]
            elif "submissionNames" in result["proteinDescription"]:
                # some entries, such as UniProtKB:A0A0G2JMH6 don't have recommendedName in proteinDescription, but follow this pattern: result->proteinDescription->submissionNames->List[0: fullName -> value].
                # there can be multiple proposed descriptions, this code accounts for them all:
                description = ""
                submissionNames = result["proteinDescription"]["submissionNames"]
                for i in range(len(submissionNames)):
                    if i == 0:
                        description = submissionNames[i]["fullName"]["value"]
                    else:
                        description += f", {submissionNames[i]['fullName']['value']}"
                # resulting description is the accumulation of all comma-delimited descriptions
            else:
                description = "ERROR: Couldn't fetch description."
                logger.warning(
                    "proteinDescription, recommendedName, fullName or value not found"
                    f" when querying for uniprot info for the id: {uniprot_id}"
                )
                logger.warning(f"result: {result}")
            (
                ensg_id,
                enst_id,
                refseq_nt_id,
            ) = self._return_mane_select_values_from_uniprot_query(result)
            return {
                "genename": name,
                "description": description,
                "ensg_id": ensg_id,
                "enst_id": enst_id,
                "refseq_nt_id": refseq_nt_id,
            }

    async def get_uniprot_info_async(
        self, uniprot_id: str, session: aiohttp.ClientSession
    ) -> dict:
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.

        Parameters:
          - (str) uniprot_id

        If the query is successful, returns the following dictionary:
            {
                "genename": name,
                "description": description,
                "ensg_id": ensg_id,
                "enst_id": enst_id,
                "refseq_nt_id": refseq_nt_id
            }

        This function automatically uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections

        Algorithm:
        This function constructs a uniprot url (https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name)
        and queries the response. A query response example: https://pastebin.com/qrWck3QG. The structure of the query response is:
        {
            "results":[
                {
                    "entryType": "UniProtKB reviewed (Swiss-Prot)",
                    "primaryAccession": "Q9NY91",
                    "organism": {
                        "scientificName": "Homo Sapiens",
                        "commonName": "Human",
                        "taxonId": 9606
                        "lineage": [...]
                    },
                    "proteinDescription": {
                        "recommendedName": {
                            "fullName": {
                                "value": "Probable glucose sensor protein SLC5A4"
                                "evidences": [...]
                            }
                        },
                        "alternativeNames": {
                            "fullName": {
                                ...
                            }
                        }
                    },
                    "genes": [
                        {
                            "genename": {
                                "value": "SLC5A4",
                                "evidences": [...]
                            },
                            "synonyms": [
                                {
                                    "value": "SAAT1"
                                },
                                {
                                    "value": "SGLT3"
                                    "evidences": [...]
                                }
                            ]
                        }
                    ],
                    "uniProtKBCrossReferences": [
                        {
                            "database": "RefSeq",
                            "id": "NP_055042.1",
                            "properties": [
                                {
                                    "key": "NucleotideSequenceId",
                                    "value": "NM_014227.2"
                                },
                                {
                                    "database": "Ensembl",
                                    "id": "ENST00000266086.6",
                                    "properties": [
                                        {
                                            "key": "ProteinId",
                                            "value": "ENSP00000266086.3"
                                        },
                                        {
                                            "key": "GeneId",
                                            "value": "ENSG00000100191.6"
                                        }
                                    ]
                                },
                                {
                                    "database": "MANE-Select",
                                    "id": "ENST00000266086.6",
                                    "properties": [
                                        {
                                            "key": "ProteinId",
                                            "value": "ENSP00000266086.3"
                                        },
                                        {
                                            "key": "RefSeqNucleotideId",
                                            "value": "NM_014227.3"
                                        },
                                        {
                                            "key": "RefSeqProteinId",
                                            "value": "NP_055042.1"
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                },
                {
                # repeat above structure for more results
                }
            ]
        }
        """
        # Extract UniProt ID if given in "database:identifier" format
        if ":" in uniprot_id:
            uniprot_id = uniprot_id.split(":")[1]

        # Attempt to cache previous function result
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_info_async.__name__}][uniprot_id={uniprot_id}]"
        previous_result = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_result is not None:
            logger.debug(
                f"Returning cached info for uniprot id {uniprot_id}: {previous_result}"
            )
            return previous_result

        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            await asyncio.sleep(self.async_request_sleep_delay)
            QUERY_RETRIES = 3  # TODO: make parameter
            i = 0
            for _ in range(QUERY_RETRIES):
                if i == (QUERY_RETRIES - 1):
                    return None
                i += 1
                try:
                    response = await session.get(url, timeout=5)
                    response_json = await response.json()
                    Cacher.store_data("url", url, response_json)
                # except(requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ServerDisconnectedError, aiohttp.ClientResponseError) as e:
                except Exception as e:
                    logger.warning(
                        f"Exception when querying info for {uniprot_id}. Exception:"
                        f" {str(e)}"
                    )
                    self.uniprot_query_exceptions.append({f"{uniprot_id}": f"{str(e)}"})
                    await asyncio.sleep(self.async_request_sleep_delay)  # sleep before retrying
                    continue

        # single query retry
        # try:
        #    response = await session.get(url, timeout=5)
        # except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError) as e:
        #    logger.warning(f"Exception when querying info for {uniprot_id}. Exception: {str(e)}")
        #    self.uniprot_query_exceptions.append({f"{uniprot_id}": f"{str(e)}"})
        #    return None

        results = response_json["results"]
        return_value = self._process_uniprot_info_query_results(results, uniprot_id)

        if return_value is not None:
            Cacher.store_data("uniprot", uniprot_data_key, return_value)

        return return_value

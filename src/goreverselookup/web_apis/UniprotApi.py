import requests
from requests.adapters import HTTPAdapter, Retry
import aiohttp
import asyncio
import json
import time
import zlib
import re
import json
import sys
from xml.etree import ElementTree

from urllib.parse import urlparse, parse_qs, urlencode

from ..util.CacheUtil import Cacher

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


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
        self.idmapping_max_retries = 30

    def get_uniprot_id(self, gene_name, taxon, get_url_only=False):
        """
        Given a gene name, returns the corresponding UniProt ID using the UniProt API.

        Parameters:
        - gene_name (str): name of the gene to search for.
        - taxon (str): NCBITaxon:xxxx or just the NCBITaxon number (xxxx)
        - TODO: retries (int): maximum number of times to retry the request in case of network errors.
        - TODO: timeout (int): timeout in seconds for the request.
        - get_url_only: only return the query url without performing the query

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections

        Returns:
        - str: UniProt ID if found, None otherwise OR the query url, if get_url_only is True
        """
        if ":" in taxon:
            taxon = taxon.split(":")[1]

        # data key is in the format [class_name][function_name][function_params]
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_id.__name__}][gene_name={gene_name}]"
        previous_uniprot_id = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_uniprot_id is not None:
            logger.debug(
                f"Cached uniprot id {previous_uniprot_id} for gene name {gene_name}"
            )
            return previous_uniprot_id

        # Define the URL to query the UniProt API
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:{taxon}&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl"

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

    def get_uniprot_info(self, uniprot_id: str, taxon_id_num: str) -> dict:
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.

        Parameters:
          - (str) uniprot_id
          - (str) taxon_id_num: NCBITaxon id number or a full NCBITaxon

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

        if ":" in taxon_id_num:
            taxon_id_num = taxon_id_num.split(":")[1]

        # Attempt to return previously cached function return value
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_info.__name__}][uniprot_id={uniprot_id},taxon={taxon_id_num}]"
        previous_info = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_info is not None:
            logger.debug(
                f"Returning cached info for uniprot id {uniprot_id}: {previous_info}"
            )
            return previous_info

        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:{taxon_id_num}&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"

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
            if result is None:
                return {}
            # name = result.get(["genes"][0]["geneName"]["value"], None)  # gene name
            # name = result.get("genes", [])[0].get("geneName", {}).get("value", None) # this creates index out of range error
            name = None
            _n = result.get("genes",None)
            if _n is not None:
                name = _n[0].get("geneName", {}).get("value", None)
            if name is None:
                return {}
            
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
        self,
        uniprot_id: str,
        session: aiohttp.ClientSession,
        organism_taxon_id_num="9606",
    ) -> dict:
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.

        Parameters:
          - (str) uniprot_id: The UniProtKB identifier
          - (int or str) organism_taxon_id_num: The number of the NCBITaxon for this uniprot id (eg. 9606 for Homo Sapiens). Also accepts full NCBITaxon:xxxx ids

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

        if ":" in organism_taxon_id_num:
            organism_taxon_id_num = organism_taxon_id_num.split(":")[1]

        # Attempt to cache previous function result
        uniprot_data_key = f"[{self.__class__.__name__}][{self.get_uniprot_info_async.__name__}][uniprot_id={uniprot_id},taxon={organism_taxon_id_num}]"
        previous_result = Cacher.get_data("uniprot", uniprot_data_key)
        if previous_result is not None:
            logger.debug(
                f"Returning cached info for uniprot id {uniprot_id}: {previous_result}"
            )
            return previous_result

        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:{organism_taxon_id_num}&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            await asyncio.sleep(self.async_request_sleep_delay)
            QUERY_RETRIES = 3
            i = 0
            for _ in range(QUERY_RETRIES):
                if i == (QUERY_RETRIES - 1):
                    return None
                i += 1
                try:
                    response = await session.get(url, timeout=5)
                    # response_json = await response.json()
                    response_content = await response.read()
                    response_json = json.loads(response_content)
                    Cacher.store_data("url", url, response_json)
                # except(requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ServerDisconnectedError, aiohttp.ClientResponseError) as e:
                except Exception as e:
                    logger.warning(
                        "Exception when querying info for UniProt ID %s: [%s] %s",
                        uniprot_id,
                        type(e).__name__,
                        str(e) if str(e) else "<no message>"
                    )
                    self.uniprot_query_exceptions.append({f"{uniprot_id}": f"{str(e)}"})
                    await asyncio.sleep(self.async_request_sleep_delay)  # sleep before retrying
                    continue

        results = response_json["results"]
        return_value = self._process_uniprot_info_query_results(results, uniprot_id)

        if return_value is not None:
            Cacher.store_data("uniprot", uniprot_data_key, return_value)

        return return_value

    def get_allowed_to_and_from_databases(self):
        """
        Functions like (UniProtApi).idmapping_batch rely on database identifiers to construct the query url to perform a batch identifier mapping.
        The specific database identifiers (db ids) that can be used as a "from_database" and "to_database" ids can be viewed using https://rest.uniprot.org/configure/idmapping/fields.

        This function returns a dictionary with two keys: 'from_dbs' and 'to_dbs', and each key holds a list of associated db ids that can be used as from_db or to_db respectively.
        """
        if hasattr(self, "idmapping_to_and_from_db_ids"):
            if (
                self.idmapping_to_and_from_db_ids is not None
                and self.idmapping_to_and_from_db_ids != {}
            ):
                return self.idmapping_to_and_from_db_ids

        url = "https://rest.uniprot.org/configure/idmapping/fields"
        r = self.s.get(url)
        r.raise_for_status()
        response = r.json()
        db_groups = response["groups"]

        # response json format:
        # {
        #    "groups": [
        #        {
        #            "groupName": "UniProt",
        #            "items": [
        #            {
        #                "displayName": "UniProtKB",
        #                "name": "UniProtKB",
        #                "from": false,
        #                "to": true,
        #                "ruleId": null,
        #                "uriLink": "https://www.uniprot.org/uniprot/%id/entry"
        #            },
        #            ...
        #            ]
        #        }
        #    "rules": ...
        # }

        # iterate database groups - currently, database groups are: "UniProt", "Sequence databases", "3D structure databases", "Protein-protein interaction databases", "Chemistry",
        # "Protein family/group databases", "PTM databases", "Genetic variation databases", "2D gel databases", "Proteomic databases", "Protocols and materials databases", "Genome annotation databases",
        # "Organism-specific databases", "Phylogenomic databases", "Enzyme and pathway databases", "Miscellaneous", "Gene expression databases" and "Family and domain databases"
        to_dbs = []
        from_dbs = []
        for db_group in db_groups:
            db_group_name = db_group.get("groupName")
            db_group_items = db_group.get("items")
            for database in db_group_items:
                db_name = database.get("name")
                from_value = bool(database.get("from"))
                to_value = bool(database.get("to"))

                if to_value == True:
                    to_dbs.append(db_name)
                if from_value == True:
                    from_dbs.append(db_name)

        self.idmapping_to_and_from_db_ids = {"from_dbs": from_dbs, "to_dbs": to_dbs}

        return self.idmapping_to_and_from_db_ids

    def idmapping_batch(
        self,
        ids: list,
        from_db: str = "UniProtKB_AC-ID",
        to_db: str = "Ensembl",
        chunk_size: int = 1000,
    ):
        """
        Performs identifier mapping from 'from_db' to 'to_db' on every identifier from 'ids',
        submitting to UniProt in sequential chunks (default size=1000) and aggregating results.

        WARNING:
        - For MGI ids (eg. MGI:1932410), UniProt expects the FULL id (MGI:1932410).
        - For ZFIN ids (eg. ZFIN:ZDB-GENE-040426-2477), UniProt expects processed ids
            (without the 'ZFIN:' prefix): ZDB-GENE-040426-2477.
        """
        # --- helpers --------------------------------------------------------------
        def check_to_and_from_db_validity(from_db, to_db):
            valid_db_ids = self.get_allowed_to_and_from_databases()
            valid_db_ids_from = valid_db_ids.get("from_dbs")
            valid_db_ids_to = valid_db_ids.get("to_dbs")

            if from_db not in valid_db_ids_from:
                raise Exception(
                    f"The value of 'from_db' {from_db} is not among allowed 'from databases' for idmapping. "
                    f"Allowed 'from databases' are: {valid_db_ids_from}"
                )
            if to_db not in valid_db_ids_to:
                raise Exception(
                    f"The value of 'to_db' {to_db} is not among allowed 'to databases' for idmapping. "
                    f"Allowed 'to databases' are: {valid_db_ids_to}"
                )
            return True

        def _chunks(seq, n):
            for i in range(0, len(seq), n):
                yield seq[i : i + n]

        def check_idmapping_results_ready(
            job_id: str,
            session: requests.Session,
            api_url: str = "https://rest.uniprot.org",
            max_retries: int = 30,
            polling_interval: float = 3.0,
        ) -> bool:
            """
            Returns True when job is FINISHED and results can be fetched.
            Returns False if we ran out of retries while still pending/running.
            Raises on FAILED/ERROR.
            """
            for i in range(1, max_retries + 1):
                r = session.get(f"{api_url}/idmapping/status/{job_id}", timeout=15)
                r.raise_for_status()
                payload = r.json()
                status = (payload.get("jobStatus") or "").upper()

                if status in ("NEW", "PENDING", "QUEUED", "RUNNING"):
                    logger.info(
                        f"idmapping for job id {job_id} in progress (retry {i}/{max_retries}), "
                        f"retrying in {polling_interval}s"
                    )
                    time.sleep(polling_interval)
                    continue

                if status in ("FINISHED", "COMPLETE"):
                    return True

                if status in ("FAILED", "ERROR"):
                    detail = payload.get("messages") or payload
                    raise Exception(
                        f"UniProt idmapping failed. jobStatus={status}; details={detail}"
                    )

                # Fallback for odd payloads that directly include results
                if "results" in payload or "failedIds" in payload:
                    return True

                logger.warning(
                    f"Unexpected status payload for job {job_id}: {payload}. "
                    f"Retrying in {polling_interval}s"
                )
                time.sleep(polling_interval)

            return False

        def get_idmapping_results_link(
            job_id: str, session: requests.Session, api_url: str = "https://rest.uniprot.org"
        ) -> str:
            url = f"{api_url}/idmapping/details/{job_id}"
            request = session.get(url, timeout=30)
            request.raise_for_status()
            j = request.json()
            return j["redirectURL"]

        def decode_results(response, file_format, compressed):
            if compressed:
                decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
                if file_format == "json":
                    return json.loads(decompressed.decode("utf-8"))
                elif file_format == "tsv":
                    return [
                        line for line in decompressed.decode("utf-8").split("\n") if line
                    ]
                elif file_format == "xlsx":
                    return [decompressed]
                elif file_format == "xml":
                    return [decompressed.decode("utf-8")]
                else:
                    return decompressed.decode("utf-8")
            elif file_format == "json":
                return response.json()
            elif file_format == "tsv":
                return [line for line in response.text.split("\n") if line]
            elif file_format == "xlsx":
                return [response.content]
            elif file_format == "xml":
                return [response.text]
            return response.text

        def get_next_link(headers):
            re_next_link = re.compile(r'<(.+)>; rel="next"')
            if "Link" in headers:
                match = re_next_link.match(headers["Link"])
                if match:
                    return match.group(1)

        def get_batch(batch_response, file_format, compressed, session):
            batch_url = get_next_link(batch_response.headers)
            while batch_url:
                batch_response = session.get(batch_url, timeout=60)
                batch_response.raise_for_status()
                yield decode_results(batch_response, file_format, compressed)
                batch_url = get_next_link(batch_response.headers)

        def print_progress_batches(batch_index, size, total):
            n_fetched = min((batch_index + 1) * size, total)
            logger.info(f"Fetched: {n_fetched} / {total}")

        def merge_xml_results(xml_results):
            merged_root = ElementTree.fromstring(xml_results[0])
            for result in xml_results[1:]:
                root = ElementTree.fromstring(result)
                for child in root.findall("{http://uniprot.org/uniprot}entry"):
                    merged_root.insert(-1, child)
            ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
            return ElementTree.tostring(
                merged_root, encoding="utf-8", xml_declaration=True
            )

        def get_xml_namespace(element):
            m = re.match(r"\{(.*)\}", element.tag)
            return m.groups()[0] if m else ""

        def combine_batches(all_results, batch_results, file_format):
            if file_format == "json":
                for key in ("results", "failedIds"):
                    if key in batch_results and batch_results[key]:
                        all_results[key] += batch_results[key]
            elif file_format == "tsv":
                return all_results + batch_results[1:]
            else:
                return all_results + batch_results
            return all_results

        def get_idmapping_results_search(url, session: requests.Session):
            parsed = urlparse(url)
            query = parse_qs(parsed.query)
            file_format = query["format"][0] if "format" in query else "json"
            if "size" in query:
                size = int(query["size"][0])
            else:
                size = 500  # page size for pagination over results
                query["size"] = size
            compressed = (
                query["compressed"][0].lower() == "true"
                if "compressed" in query
                else False
            )
            parsed = parsed._replace(query=urlencode(query, doseq=True))
            url = parsed.geturl()
            request = session.get(url, timeout=60)
            try:
                results = decode_results(request, file_format, compressed)
            except json.JSONDecodeError as e:
                logger.error(f"JSONDecode error during uniprot batch request: {e}")
                return None

            total_header = request.headers.get("x-total-results")
            total = int(total_header) if total_header and total_header.isdigit() else 0
            if total:
                print_progress_batches(0, size, total)

            for i, batch in enumerate(get_batch(request, file_format, compressed, session), 1):
                results = combine_batches(results, batch, file_format)
                if total:
                    print_progress_batches(i, size, total)

            if file_format == "xml":
                return merge_xml_results(results)
            return results

        # --- validation & preprocessing ------------------------------------------
        logger.info(f"Starting UniProt idmapping batch: from {from_db} to {to_db} for {len(ids)} ids.")
        
        check_to_and_from_db_validity(from_db=from_db, to_db=to_db)

        processed_ids = []
        for _id in ids:
            if "MGI:" in _id:
                processed_ids.append(_id)  # keep full MGI:xxxx
                continue
            processed_ids.append(_id.split(":")[1]) if ":" in _id else processed_ids.append(_id)

        new_ids = processed_ids

        polling_interval = 3.0  # seconds between polls
        api_url = "https://rest.uniprot.org"

        retries = Retry(
            total=5,
            backoff_factor=0.3,
            status_forcelist=[429, 500, 502, 503, 504],
            respect_retry_after_header=True,
        )
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))

        # --- aggregate over chunks -----------------------------------------------
        all_results = {"results": [], "failedIds": []}

        for chunk_idx, chunk in enumerate(_chunks(new_ids, chunk_size), start=1):
            logger.info(
                f"Submitting UniProt idmapping chunk {chunk_idx}: size={len(chunk)} "
                f"(from={from_db} to={to_db})"
            )
            try:
                # use the SAME retrying session for submit/poll/fetch
                req = session.post(
                    f"{api_url}/idmapping/run",
                    data={"from": from_db, "to": to_db, "ids": ",".join(chunk)},
                    timeout=60,
                )
                logger.info(
                    f"Attempting batch idmapping from {from_db} to {to_db}. Url = {req.url}"
                )
                req.raise_for_status()
                job_id = req.json()["jobId"]
            except Exception as e:
                logger.error(
                    f"Error submitting UniProt idmapping for chunk {chunk_idx}: {e}"
                )
                # If a submission fails outright, mark the whole chunk as failed
                all_results["failedIds"].extend(chunk)
                continue

            try:
                ready = check_idmapping_results_ready(
                    job_id,
                    session,
                    api_url=api_url,
                    max_retries=self.idmapping_max_retries,
                    polling_interval=polling_interval,
                )
            except Exception as e:
                logger.error(
                    f"UniProt idmapping reported failure for chunk {chunk_idx} (job {job_id}): {e}"
                )
                all_results["failedIds"].extend(chunk)
                continue

            if not ready:
                logger.error(
                    f"UniProt idmapping did not finish within {self.idmapping_max_retries} polls "
                    f"for chunk {chunk_idx} (job {job_id})."
                )
                all_results["failedIds"].extend(chunk)
                continue

            try:
                link = get_idmapping_results_link(job_id, session, api_url=api_url)
                results = get_idmapping_results_search(link, session)
            except requests.exceptions.SSLError as e:
                logger.warning(
                    f"SSL exception encountered when fetching results for chunk {chunk_idx} "
                    f"(job {job_id}). Error: {e}"
                )
                results = None
            except Exception as e:
                logger.error(
                    f"Error fetching/decoding results for chunk {chunk_idx} (job {job_id}): {e}"
                )
                results = None

            if results is None:
                # if we couldn't fetch results, count the whole chunk as failed
                all_results["failedIds"].extend(chunk)
                continue

            successful = results.get("results", []) if isinstance(results, dict) else []
            failed_ids = results.get("failedIds", []) if isinstance(results, dict) else []

            all_results["results"].extend(successful)
            all_results["failedIds"].extend(failed_ids)

            logger.info(
                f"Chunk {chunk_idx} done: +{len(successful)} mapped, +{len(failed_ids)} failed"
            )

        return all_results

    def idmapping_ensembl_batch(self, uniprot_ids: list):
        """
        Maps a list of uniprot ids to a list of ensembl ids.
        """
        # possible from and to parameters: https://jsonformatter.org/f5632a
        from_db = "UniProtKB_AC-ID"
        to_db = "Ensembl"
        return self.idmapping_batch(ids=uniprot_ids, from_db=from_db, to_db=to_db)

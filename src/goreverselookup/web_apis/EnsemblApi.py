import requests
from requests.adapters import HTTPAdapter, Retry
import aiohttp
import asyncio

from ..util.CacheUtil import Cacher

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class EnsemblApi:
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
        self.ortholog_query_exceptions = (
            []
        )  # the list of exceptions during the ortholog query
        self.async_request_sleep_delay = 0.5

    def get_human_ortholog(self, id: str):
        """
        Given a source ID, detect organism and returns the corresponding human ortholog using the Ensembl API.

        Parameters:
          - (str) id

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections
        """
        ensembl_data_key = (
            f"[{self.__class__.__name__}][{self.get_human_ortholog.__name__}][id={id}]"
        )
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result

        full_id = id
        if "ZFIN" in id:
            species = "zebrafish"
            id_url = id.split(":")[1]
        elif "Xenbase" in id:
            species = "xenopus_tropicalis"
            id_url = id.split(":")[1]
        elif "MGI" in id:
            species = "mouse"
            id_url = id
        elif "RGD" in id:
            species = "rat"
            id_url = id.split(":")[1]
        else:
            logger.info(f"No predefined organism found for {id}")
            return None

        url = f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if (
            previous_response is not None
        ):  # "error" check is a bugfix for this response: {'error': 'No valid lookup found for symbol Oxct2a'}
            response_json = previous_response
        else:
            try:
                response = self.s.get(
                    url, headers={"Content-Type": "application/json"}, timeout=5
                )
                response.raise_for_status()
                response_json = response.json()["data"][0]["homologies"]
                # ConnectionCacher.store_url(url, response=response_json)
                Cacher.store_data("url", url, response_json)
            except requests.exceptions.RequestException:
                return None

        if response_json == []:
            return None

        max_perc_id = 0.0
        ortholog = ""
        for ortholog_dict in response_json:
            if ortholog_dict["target"]["species"] == "homo_sapiens":
                current_perc_id = ortholog_dict["target"]["perc_id"]
                if current_perc_id > max_perc_id:
                    ortholog = ortholog_dict["target"]["id"]
                    max_perc_id = current_perc_id

        # above code also accounts for the homo sapiens species
        # best_ortholog_dict = max(response_json, key=lambda x: int(x["target"]["perc_id"]))
        # ortholog = best_ortholog_dict["target"].get("id")

        Cacher.store_data("ensembl", ensembl_data_key, ortholog)
        logger.info(f"Received ortholog for id {full_id} -> {ortholog}")
        return ortholog

    async def get_human_ortholog_async(self, id, session: aiohttp.ClientSession):
        """
        Given a source ID, detect organism and returns the corresponding human ortholog using the Ensembl API.
        Example source IDs are: UniProtKB:P21709, RGD:6494870, ZFIN:ZDB-GENE-040426-1432, Xenbase:XB-GENE-479318 and MGI:95537

        Parameters:
          - (str) id

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections
        """
        ensembl_data_key = f"[{self.__class__.__name__}][{self.get_human_ortholog_async.__name__}][id={id}]"
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result

        if "ZFIN" in id:
            species = "zebrafish"
            id_url = id.split(":")[1]
        elif "Xenbase" in id:
            species = "xenopus_tropicalis"
            id_url = id.split(":")[1]
        elif "MGI" in id:
            species = "mouse"
            id_url = id
        elif "RGD" in id:
            species = "rat"
            id_url = id.split(":")[1]
        else:
            logger.info(f"No predefined organism found for {id}")
            return None
        url = f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            try:
                response = await session.get(
                    url, headers={"Content-Type": "application/json"}, timeout=10
                )
                # response.raise_for_status()
                response_json = await response.json()
                Cacher.store_data("url", url, response_json)
                await asyncio.sleep(self.async_request_sleep_delay)
            # except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ClientResponseError) as e:
            except Exception as e:
                logger.warning(
                    f"Exception for {id_url} for request:"
                    f" https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none."
                    f" Exception: {str(e)}"
                )
                self.ortholog_query_exceptions.append({f"{id}": f"{str(e)}"})
                return None

        # TODO: implement this safety check, server may send text only, which causes error (content_type == "text/plain")
        # if response.content_type == "application/json":
        #    response_json = await response.json()

        if response_json == [] or "error" in response_json:
            return None
        elif response_json != [] and "error" not in response_json:
            response_json = response_json["data"][0]["homologies"]
            if response_json == []:  # if there are no homologies, return None
                return None

            max_perc_id = 0.0
            ortholog = ""
            for ortholog_dict in response_json:
                if ortholog_dict["target"]["species"] == "homo_sapiens":
                    current_perc_id = ortholog_dict["target"]["perc_id"]
                    if current_perc_id > max_perc_id:
                        ortholog = ortholog_dict["target"]["id"]
                        max_perc_id = current_perc_id

            # Above code is better, because it accounts for if the "species" in the response is "homo_sapiens"
            # best_ortholog_dict = max(response_json, key=lambda x: int(x["target"]["perc_id"]))
            # ortholog = best_ortholog_dict["target"].get("id")

            Cacher.store_data("ensembl", ensembl_data_key, ortholog)
            logger.info(f"Received ortholog for id {id} -> {ortholog}")
            return ortholog

    def get_sequence(self, ensembl_id, sequence_type="cdna"):
        """
        Given an Ensembl ID, returns the corresponding nucleotide sequence using the Ensembl API.
        type can be genomic,cds,cdna,protein
        """
        url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?object_type=transcript;type={sequence_type}"
        try:
            response = self.s.get(
                url, headers={"Content-Type": "text/plain"}, timeout=5
            )
            response.raise_for_status()
        except requests.exceptions.RequestException:
            logger.warning(f"Failed to fetch Ensembl sequence for {ensembl_id}")
            return None
        sequence = response.text
        logger.info(f"Received sequence for id {ensembl_id}.")
        return sequence

    def get_info(self, id: str) -> dict:
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol

        Returns:
            dict: Information about the gene
        """
        if (
            "Error" in id
        ):  # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(
                f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene"
                " does not have a human ortholog and you are safe to ignore it."
            )
            return {}

        ensembl_data_key = (
            f"[{self.__class__.__name__}][{self.get_info.__name__}][id={id}]"
        )
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached info for id {id}: {previous_result}")
            return previous_result

        species_mapping = {
            "ZFIN": "zebrafish/",
            "Xenbase": "xenopus_tropicalis/",
            "MGI": "mouse/MGI:",
            "RGD": "rat/",
            "UniProtKB": "human/",
        }

        # Check if the ID is an Ensembl ID or symbol
        if id.startswith("ENS"):
            endpoint = f"id/{id}"
        else:
            # TODO: Check if this is ever even queried. I am trying to see if any such urls are stored in connection_cache.json, but I see none.
            #
            # One of the following links is used:
            #   rest.ensembl.org/lookup/symbol/zebrafish/{ZFIN_ID}
            #   rest.ensembl.org/lookup/symbol/xenopus_tropicalis/{XENBASE_ID}
            #   rest.ensembl.org/lookup/symbol/mouse/{MGI_ID}
            #   rest.ensembl.org/lookup/symbol/rat/{RGD_ID}
            #   rest.ensembl.org/lookup/symbol/human/{UNIPROT_ID}

            prefix, id_ = id.split(":") if ":" in id else (None, id)
            species = species_mapping.get(
                prefix, "human/"
            )  # defaults to human if not prefix "xxx:"
            endpoint = f"symbol/{species}{id_}"

        url = f"https://rest.ensembl.org/lookup/{endpoint}?mane=1;expand=1"

        try:
            # Check if the url is cached
            # previous_response = ConnectionCacher.get_url_response(url)
            previous_response = Cacher.get_data("url", url)
            if previous_response is not None:
                response_json = previous_response
            else:
                response = self.s.get(
                    url,
                    headers={"Content-Type": "application/json"},
                    timeout=5,
                )
                response.raise_for_status()
                response_json = response.json()
                # ConnectionCacher.store_url(url, response_json)
                Cacher.store_data("url", url, response_json)
        except requests.exceptions.RequestException:
            # If the request fails, try the xrefs URL instead
            try:
                if (
                    "ENS" not in id
                ):  # id is not an ensembl id, attempt to find a cross-reference
                    url = f"https://rest.ensembl.org/xrefs/{endpoint}?"
                    # Check if the url is cached
                    previous_response = Cacher.get_data("url", url)
                    if previous_response is not None:
                        response_json = previous_response
                    else:
                        response = self.s.get(
                            url,
                            headers={"Content-Type": "application/json"},
                            timeout=5,
                        )
                        response.raise_for_status()
                        response_json = response.json()
                        Cacher.store_data("url", url, response_json)
                    # Use the first ENS ID in the xrefs response to make a new lookup request
                    ensembl_id = next(
                        (xref["id"] for xref in response_json if "ENS" in xref["id"]),
                        None,
                    )
                else:
                    # id is an ensembl id
                    ensembl_id = id

                if ensembl_id:
                    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?mane=1;expand=1"
                    # Check if the url is cached
                    # previous_response = ConnectionCacher.get_url_response(url)
                    previous_response = Cacher.get_data("url", url)
                    if previous_response is not None:
                        response_json = previous_response
                    else:
                        response = self.s.get(
                            url,
                            headers={"Content-Type": "application/json"},
                            timeout=5,
                        )
                        response.raise_for_status()
                        response_json = response.json()
                        # ConnectionCacher.store_url(url, response_json)
                        Cacher.store_data("url", url, response_json)
                else:
                    raise Exception("no ensembl id returned")
            except Exception:
                logger.warning(f"Failed to fetch Ensembl info for {id}.")
                return {}

        # Extract gene information from API response
        ensg_id = response_json.get("id")
        name = response_json.get("display_name")
        description = response_json.get("description", "").split(" [")[0]

        canonical_transcript_id = next(
            (
                entry.get("id")
                for entry in response_json["Transcript"]
                if entry.get("is_canonical")
            ),
            None,
        )
        mane_transcripts = [d for d in response_json["Transcript"] if d.get("MANE")]
        if len(mane_transcripts) == 0:
            ensembl_transcript_id = canonical_transcript_id
            refseq_id = None
        elif len(mane_transcripts) == 1:
            ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get("id")
            refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
        else:
            selected_entry = next(
                (entry for entry in mane_transcripts if entry.get("is_canonical")), None
            )
            if not selected_entry:
                ensembl_transcript_id = selected_entry["MANE"][0].get("id")
                refseq_id = selected_entry["MANE"][0].get("refseq_match")
            else:
                ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get(
                    "id"
                )  # select the first canonical transcript with MANE
                refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
                logger.warning(f"Found non-canonical MANE transcript for {id}")

        if ensembl_transcript_id:
            try:
                url = f"https://rest.ensembl.org/xrefs/id/{ensembl_transcript_id}?all_levels=1;external_db=UniProt%"
                # previous_response = ConnectionCacher.get_url_response(url)
                previous_response = Cacher.get_data("url", url)
                if previous_response is not None:
                    response_json = previous_response
                else:
                    response = self.s.get(
                        url,
                        headers={"Content-Type": "application/json"},
                        timeout=5,
                    )
                    response.raise_for_status()
                    response_json = response.json()
                    # ConnectionCacher.store_url(url, response_json)
                    Cacher.store_data("url", url, response_json)
            except requests.exceptions.RequestException:
                pass
            uniprot_id = ""
            # bugfix: attribute error, because some 'entry' objects in loop were read as strings
            # uniprot_id = next((entry.get("primary_id") for entry in response_json if entry.get("dbname") =="Uniprot/SWISSPROT"), None)
            for entry in response_json:
                if isinstance(entry, dict):
                    if entry.get("dbname") == "Uniprot/SWISSPROT":
                        uniprot_id = entry.get("primary_id")

        logger.debug(f"Received info data for id {id}.")
        return_value = {
            "ensg_id": ensg_id,
            "genename": name,
            "description": description,
            "enst_id": ensembl_transcript_id,
            "refseq_nt_id": refseq_id,
            "uniprot_id": uniprot_id,
        }
        Cacher.store_data("ensembl", ensembl_data_key, return_value)
        return return_value

    async def get_info_async(self, id: str, session: aiohttp.ClientSession):
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol

        Returns:
            dict: Information about the gene
        """
        if ("Error" in id):  # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(
                f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene"
                " does not have a human ortholog and you are safe to ignore it."
            )
            return {}

        ensembl_data_key = (
            f"[{self.__class__.__name__}][{self.get_info_async.__name__}][id={id}]"
        )
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result

        species_mapping = {
            "ZFIN": "zebrafish/",
            "Xenbase": "xenopus_tropicalis/",
            "MGI": "mouse/MGI:",
            "RGD": "rat/",
            "UniProtKB": "human/",
        }

        # Check if the ID is an Ensembl ID or symbol
        if id.startswith("ENS"):
            endpoint = f"id/{id}"
        else:
            prefix, id_ = id.split(":") if ":" in id else (None, id)
            species = species_mapping.get(
                prefix, "human/"
            )  # defaults to human if not prefix "xxx:"
            endpoint = f"symbol/{species}{id_}"

        try:
            # TODO: 19.08.2023: the below link doesn't work for any other {species} in endpoint other than human. Ie.
            # zebrafish, xenbase, mgi, rgd don't work !!!
            url = f"https://rest.ensembl.org/lookup/{endpoint}?mane=1;expand=1"
            # previous_response = ConnectionCacher.get_url_response(url)
            previous_response = Cacher.get_data("url", url)
            if previous_response is not None:
                response_json = previous_response
            else:
                response = await session.get(
                    url, headers={"Content-Type": "application/json"}, timeout=5
                )
                # response.raise_for_status()
                response_json = await response.json()
                # ConnectionCacher.store_url(url, response_json)
                Cacher.store_data("url", url, response_json)
                await asyncio.sleep(self.async_request_sleep_delay)
        except (
            requests.exceptions.RequestException,
            TimeoutError,
            asyncio.CancelledError,
            asyncio.exceptions.TimeoutError,
            aiohttp.ClientResponseError,
            aiohttp.ClientOSError,
            aiohttp.ClientPayloadError,
        ):
            # If the request fails, try the xrefs URL instead
            try:
                # TODO: 19.08.2023: the below link doesn't work for any other {species} in endpoint other than human. Ie.
                # zebrafish, xenbase, mgi, rgd don't work !!!
                # The xrefs link is supposed to work for parameter 'id's which are not ENSG
                ensembl_id = ""
                if "ENS" not in id:
                    # parameter id is not ensembl, attempt to find ensembl id
                    url = (  # cross references
                        f"https://rest.ensembl.org/xrefs/{endpoint}?"
                    )
                    previous_response = Cacher.get_data("url", url)
                    if previous_response is not None:
                        response_json = previous_response
                    else:
                        response = await session.get(
                            url, headers={"Content-Type": "application/json"}, timeout=5
                        )
                        # response.raise_for_status()
                        response_json = await response.json()
                        Cacher.store_data("url", url, response_json)
                        await asyncio.sleep(self.async_request_sleep_delay)
                    # Use the first ENS ID in the xrefs response to make a new lookup request
                    ensembl_id = next(
                        (xref["id"] for xref in response_json if "ENS" in xref["id"]),
                        None,
                    )
                else:
                    # ensembl id is already supplied in the parameter id
                    ensembl_id = id

                if ensembl_id:
                    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?mane=1;expand=1"
                    previous_response = Cacher.get_data("url", url)
                    if previous_response is not None:
                        response_json = previous_response
                    else:
                        response = await session.get(
                            url, headers={"Content-Type": "application/json"}, timeout=5
                        )
                        # response.raise_for_status()
                        response_json = await response.json()
                        Cacher.store_data("url", url, response_json)
                        await asyncio.sleep(self.async_request_sleep_delay)
                else:
                    raise Exception("no ensembl id returned")
            except Exception as e:
                logger.warning(f"Failed to fetch Ensembl info for {id}. Error = {e}")
                return {}

        if "error" in response_json or response_json is None:
            logger.warning(f"Failed to fetch Ensembl info for {id}.")
            return {}

        # Extract gene information from API response
        ensg_id = response_json.get("id")
        name = response_json.get("display_name")
        description = response_json.get("description", "").split(" [")[0]

        canonical_transcript_id = next(
            (
                entry.get("id")
                for entry in response_json["Transcript"]
                if entry.get("is_canonical")
            ),
            None,
        )
        mane_transcripts = [d for d in response_json["Transcript"] if d.get("MANE")]
        if len(mane_transcripts) == 0:
            ensembl_transcript_id = canonical_transcript_id
            refseq_id = None
        elif len(mane_transcripts) == 1:
            ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get("id")
            refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
        else:
            selected_entry = next(
                (entry for entry in mane_transcripts if entry.get("is_canonical")), None
            )
            if not selected_entry:
                ensembl_transcript_id = selected_entry["MANE"][0].get("id")
                refseq_id = selected_entry["MANE"][0].get("refseq_match")
            else:
                ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get(
                    "id"
                )  # select the first canonical transcript with MANE
                refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
                logger.warning(f"Found non-canonical MANE transcript for {id}")

        if ensembl_transcript_id:
            try:
                url = f"https://rest.ensembl.org/xrefs/id/{ensembl_transcript_id}?all_levels=1;external_db=UniProt%"
                previous_response = Cacher.get_data("url", url)
                if previous_response is not None:
                    response_json = previous_response
                else:
                    response = await session.get(
                        url, headers={"Content-Type": "application/json"}, timeout=5
                    )
                    response.raise_for_status()  # TODO: solve Too Many Requests error (429) -> aiohttp.client_exceptions.ClientResponseError: 429, message='Too Many Requests', url=URL('https://rest.ensembl.org/xrefs/id/ENST00000301012?all_levels=1;external_db=UniProt%25')
                    response_json = await response.json()
                    Cacher.store_data("url", url, response_json)
                    await asyncio.sleep(self.async_request_sleep_delay)
            # except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ClientResponseError):
            #    pass
            except Exception as e:
                logger.warning(f"Exception: {e}")
                if "Too Many Requests" in (f"{e}"):
                    logger.info(f"Too many requests detected. Sleeping for 5 seconds.")
                    await asyncio.sleep(5)
                pass

            uniprot_id = ""
            # bugfix: attribute error, because some 'entry' objects in loop were read as strings
            # uniprot_id = next((entry.get("primary_id") for entry in response_json if entry.get("dbname") =="Uniprot/SWISSPROT"), None)
            for entry in response_json:
                if isinstance(entry, dict):
                    if entry.get("dbname") == "Uniprot/SWISSPROT":
                        uniprot_id = entry.get("primary_id")

        logger.debug(f"Received info data for id {id}.")
        return_value = {
            "ensg_id": ensg_id,
            "genename": name,
            "description": description,
            "enst_id": ensembl_transcript_id,
            "refseq_nt_id": refseq_id,
            "uniprot_id": uniprot_id,
        }
        Cacher.store_data("ensembl", ensembl_data_key, return_value)
        return return_value

import requests
from requests.adapters import HTTPAdapter, Retry
import aiohttp
import asyncio
import json

from ..util.CacheUtil import Cacher
from ..util.ApiUtil import EnsemblUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger

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

        TODO: update to full taxon modularity
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
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:  # "error" check is a bugfix for this response: {'error': 'No valid lookup found for symbol Oxct2a'}
            response_json = previous_response
        else:
            try:
                response = self.s.get(
                    url, 
                    headers={"Content-Type": "application/json"}, 
                    timeout=5
                )
                response.raise_for_status()
                response_json = response.json()["data"][0]["homologies"]
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
    
    def batch_ensembl_lookup(self, ids:list, max_ids_per_request:int = 1000):
        """
        Performs a batch Ensembl lookup/id/<id> query for all input Ensembl ids.

        Params:
          - (list) ids: input Ensembl ids
          - (int) max_ids_per_request: this is a value specified by Ensembl and do not change it
        """
        num_ids = len(ids)
        logger.info(f"Received {num_ids} ids for batch Ens lookup.")
        sublists = []
        if num_ids > max_ids_per_request:
            for i in range(0, num_ids, max_ids_per_request):
                end = i+max_ids_per_request if (i+max_ids_per_request)<num_ids else num_ids # define the upper boundary
                start = i
                sublist = ids[start:end]
                sublists.append(sublist)
        else:
            sublists = [ids]
        
        processed_count = 0
        responses = []
        for sublist in sublists:
            processed_count += len(sublist)
            logger.info(f"Processing {processed_count}/{num_ids} ids. Please, be patient.")

            server = "https://rest.ensembl.org"
            ext = "/lookup/id"
            headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
            
            data = json.dumps({"ids": sublist})
            data = str(data)

            # data = f"{{ \"ids\" : {sublist} }}"
            
            #data = {}
            #data["ids"] = ids
            #data = str(data) # example: '{ "ids" : ["ENSG00000157764", "ENSG00000248378" ] }'
            r = requests.post(server+ext, headers=headers, data=data)
            if not r.ok:
                r.raise_for_status()
                #sys.exit()

            response_json = r.json()
            responses.append(response_json)
        
        merged_dict = {}
        for response in responses:
            for key, value in response.items():
                if key in merged_dict:
                    merged_dict[key].update(value)
                else:
                    merged_dict[key] = value
        
        return merged_dict

    async def get_human_ortholog_async(self, id, session: aiohttp.ClientSession, taxon=""):
        """
        Given a source ID, detect organism and returns the corresponding human ortholog using the Ensembl API.
        Example source IDs are: UniProtKB:P21709, RGD:6494870, ZFIN:ZDB-GENE-040426-1432, Xenbase:XB-GENE-479318 and MGI:95537.

        Taxon can be specified, since UniProtKB ids may not be of human origin.

        Parameters:
          - (str) id
          - (aiohttp.ClientSession) session
          - (str) taxon: A full NCBITaxon (eg. NCBITaxon:9606).
                         Currently, only the following taxa are supported: NCBITaxon:9606 (Homo Sapiens), NCBITaxon:7955 (Danio rerio),  NCBITaxon:10116 (Rattus norvegicus), NCBITaxon:10090 (Mus musculus), NCBITaxon:8353 (Xenopus)

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections
        """
        # TODO: implement modular taxon numbers!!
        def taxon_to_label(taxon:str):
            """
            Converts a NCBI taxon to a suitable ensembl label
            """
            match taxon:
                case "NCBITaxon:9606":
                    return "homo_sapiens"
                case "NCBITaxon:7955":
                    return "danio_rerio"
                case "NCBITaxon:10116":
                    return "rattus_norvegicus"
                case "NCBITaxon:10090":
                    return "mus_musculus"
                case "NCBITaxon:8353":
                    return "xenopus"

        ensembl_data_key = f"[{self.__class__.__name__}][{self.get_human_ortholog_async.__name__}][id={id}]"
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result

        id_url = None
        if "UniProtKB" in id:
            id_url = id.split(":")[1]
        elif "ZFIN" in id:
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
            # logger.info(f"No predefined organism found for {id}")
            # return None
            pass

        if id_url == None:
            # attempt one final split
            id_url = id.split(":")[1]
        
        if taxon != "" and taxon is not None:
            species = taxon_to_label(taxon=taxon)

        # TODO: SWITCH TARGET SPECIES HERE !!!
        # this link is important to query between 3rd party databases (zfin, mgi, rgd, xenbase, ...) and the target organism
        url = f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            try:
                response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=10)
                # response.raise_for_status()
                # response_json = await response.json()
                response_content = await response.read()
                response_json = json.loads(response_content)
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
            try:
                response_json = response_json["data"][0]["homologies"]
            except (KeyError, IndexError):
                logger.warning(f"Key error or list index out of range when parsing ['data'][0]['homologies] for id {id}.")
                logger.debug(f"response json = {response_json}")
                return None
            if response_json == []:  # if there are no homologies, return None
                return None

            max_perc_id = 0.0
            ortholog = ""

            # debug check if any zfin,mgi,rgd,xenbase actually make it up to here
            if species in ["zebrafish", "mouse", "rat", "xenopus_tropicalis", "danio_rerio", "mus_musculus", "rattus_norvegicus", "xenopus"]:
                pass # TODO: remove this
            
            # debug check if any non-homo-sapiens uniprotkb genes get processed
            if species in ["zebrafish", "mouse", "rat", "xenopus_tropicalis", "danio_rerio", "mus_musculus", "rattus_norvegicus", "xenopus"] and "UniProtKB" in id:
                pass # TODO: remove this

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
    
    async def get_ortholog_async(self, id, session: aiohttp.ClientSession, source_taxon:str, target_taxon:str):
        """
        Given a source ID, detect organism and returns the corresponding ortholog (for the target taxon).
        Example source IDs are: RGD:6494870, ZFIN:ZDB-GENE-040426-1432, Xenbase:XB-GENE-479318 and MGI:95537.

        Parameters:
          - (str) id: A 3rd party id (such as ZFIN:ZDB-GENE-040425-1432) or an Ensembl id (eg. ENSGxxxx)
          - (aiohttp.ClientSession) session
          - (str) taxon: A full NCBITaxon (eg. NCBITaxon:9606).

        Warning: Using UniProtKB ids always fails. See https://pastebin.com/ffMyJTBK. 

        This function uses request caching. It will use previously saved url request responses instead of performing new (the same as before) connections
        """
        ensembl_data_key = f"[{self.__class__.__name__}][{self.get_human_ortholog_async.__name__}][id={id}]"
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result
        
        id_url = None
        if "MGI" in id: # MGI ids are expected in MGI:XXXX format
            id_url = id
        else:
            # UniProtKB, ZFIN, Xenbase, RGD are expected without the starting database identifier
            if ":" in id:
                id_url = id.split(":")[1]
        
        if source_taxon is None or target_taxon is None:
            return
        source_species_label = EnsemblUtil.taxon_to_ensembl_label(source_taxon)
        target_species_label = EnsemblUtil.taxon_to_ensembl_label(target_taxon)

        # this link doesn't work for UniProtKB ids (see https://pastebin.com/ffMyJTBK)
        # works for ZDB-GENE-... ids and for Ensembl ids (example: ZDB-GENE-990708-7, OTPB_DANRE, UniProtKB:Q6DGH9, ENSDARG00000058379)
        url = f"https://rest.ensembl.org/homology/symbol/{source_species_label}/{id_url}?target_species={target_species_label};type=orthologues;sequence=none"

        # Check if the url is cached
        # previous_response = ConnectionCacher.get_url_response(url)
        previous_response = Cacher.get_data("url", url)
        if previous_response is not None:
            response_json = previous_response
        else:
            try:
                response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=15)
                # response.raise_for_status()
                # response_json = await response.json()
                response_content = await response.read()
                response_json = json.loads(response_content)
                Cacher.store_data("url", url, response_json)
                await asyncio.sleep(self.async_request_sleep_delay)
            # except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ClientResponseError) as e:
            except Exception as e:
                logger.warning(
                    f"Exception for {id_url} for request:"
                    f"{url}"
                    f" Exception: {str(e)}"
                )
                self.ortholog_query_exceptions.append({f"{id}": f"{str(e)}"})
                return None

        if response_json == [] or "error" in response_json:
            return None
        elif response_json != [] and "error" not in response_json:
            try:
                response_json = response_json["data"][0]["homologies"]
            except (KeyError, IndexError):
                logger.warning(f"Key error or list index out of range when parsing ['data'][0]['homologies] for id {id}.")
                logger.debug(f"response json = {response_json}")
                return None
            if response_json == []:  # if there are no homologies, return None
                return None

            max_perc_id = 0.0
            ortholog = ""

            for ortholog_dict in response_json:
                if ortholog_dict["target"]["species"] == target_species_label:
                    current_perc_id = ortholog_dict["target"]["perc_id"]
                    if current_perc_id > max_perc_id:
                        ortholog = ortholog_dict["target"]["id"]
                        max_perc_id = current_perc_id

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

    def get_info(self, id: str, taxon:str) -> dict:
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol
            taxon (str): NCBI taxon for the id

        Returns:
            dict: Information about the gene
        """
        if "Error" in id:  # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene does not have a human ortholog and you are safe to ignore it.")
            return {}
        
        if isinstance(id, list):
            id = id[0]
        
        # some ids can have a ".", such as ENSDARG0001245.1 -> remove what is after the .
        if "." in id:
            id = id.split(".")[0]

        ensembl_data_key = (f"[{self.__class__.__name__}][{self.get_info.__name__}][id={id},taxon={taxon}]")
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached info for id {id}: {previous_result}")
            return previous_result
        
        if ":" in id:
            # split all ids besides MGI, since MGI ids are expected in the MGI:xxxx format
            # other ids eg. ZFIN:ZDB-GENE-xxx are expected in ZDB-GENE-xxx format
            if "MGI" not in id:
                id = id.split(":")[1]
        
        if id.startswith("ENS"):
            endpoint = f"id/{id}"
        else:
            species = EnsemblUtil.taxon_to_ensembl_label(taxon_num=taxon)
            endpoint = f"smybol/{species}/{id}"

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
                if "ENS" not in id:  # id is not an ensembl id, attempt to find a cross-reference
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

    async def get_info_async(self, id: str, taxon:str, session: aiohttp.ClientSession, request_timeout = 15):
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol
            taxon (str): The taxon for this id (a NCBITaxon or a taxon nummber identifier)

        Returns:
            dict: Information about the gene
        
        WARNING: Leave request timeout at 15 or greater. When request timeout is 5, many requests fail.
        """
        if ("Error" in id):  # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(
                f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene"
                " does not have a human ortholog and you are safe to ignore it."
            )
            return {}
        
        if isinstance(id, list):
            id = id[0]
        
        # some ids can have a ".", such as ENSDARG0001245.1 -> remove what is after the .
        if "." in id:
            id = id.split(".")[0]

        ensembl_data_key = (f"[{self.__class__.__name__}][{self.get_info_async.__name__}][id={id},taxon={taxon}]")
        previous_result = Cacher.get_data("ensembl", ensembl_data_key)
        if previous_result is not None:
            logger.debug(f"Returning cached ortholog for id {id}: {previous_result}")
            return previous_result
        
        
        if ":" in id:
            # split all ids besides MGI, since MGI ids are expected in the MGI:xxxx format
            # other ids eg. ZFIN:ZDB-GENE-xxx are expected in ZDB-GENE-xxx format
            if "MGI" not in id:
                id = id.split(":")[1]
        
        if id.startswith("ENS"):
            endpoint = f"id/{id}"
        else:
            species = EnsemblUtil.taxon_to_ensembl_label(taxon_num=taxon)
            endpoint = f"smybol/{species}/{id}"

        try:
            # TODO: 19.08.2023: the below link doesn't work for any other {species} in endpoint other than human. Ie.
            # zebrafish, xenbase, mgi, rgd don't work !!!
            url = f"https://rest.ensembl.org/lookup/{endpoint}?mane=1;expand=1"
            previous_response = Cacher.get_data("url", url)
            if previous_response is not None:
                response_json = previous_response
            else:
                # headers = {"content-type": "application/json"}
                # full_url = f"{url}?{'&'.join([f'{key}={value}' for key, value in headers.items()])}"
                # logger.debug(f"Attempting url: {full_url}")
                response = await session.get(url, headers={"content-type": "application/json"}, timeout=request_timeout)
                response.raise_for_status()
                # response_json = await response.json()
                response_content = await response.read()
                response_json = json.loads(response_content)
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
            asyncio.exceptions.TimeoutError
        ) as e:
            logger.debug(f"Exception {e} with url {url}. Continuing info query with alternate urls.")
            # If the request fails, try the xrefs URL instead
            try:
                # TODO: 19.08.2023: the below link doesn't work for any other {species} in endpoint other than human. Ie.
                # zebrafish, xenbase, mgi, rgd don't work !!!
                # The xrefs link is supposed to work for parameter 'id's which are not ENSG
                ensembl_id = ""
                if "ENS" not in id:
                    # parameter id is not ensembl, attempt to find ensembl id cross references
                    url = f"https://rest.ensembl.org/xrefs/{endpoint}?"
                    previous_response = Cacher.get_data("url", url)
                    if previous_response is not None:
                        response_json = previous_response
                    else:
                        response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=request_timeout)
                        # response.raise_for_status()
                        # response_json = await response.json()
                        response_content = await response.read()
                        response_json = json.loads(response_content)
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
                        response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=request_timeout)
                        # response.raise_for_status()
                        # response_json = await response.json()
                        response_content = await response.read()
                        response_json = json.loads(response_content)
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
                    response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=5)
                    response.raise_for_status()  # TODO: solve Too Many Requests error (429) -> aiohttp.client_exceptions.ClientResponseError: 429, message='Too Many Requests', url=URL('https://rest.ensembl.org/xrefs/id/ENST00000301012?all_levels=1;external_db=UniProt%25')
                    # response_json = await response.json()
                    response_content = await response.read()
                    response_json = json.loads(response_content)
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

        logger.debug(f"Received info data for id {id}: ensg_id = {ensg_id}, genename = {name}, enst_id = {ensembl_transcript_id}, uniprot_id = {uniprot_id}")
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

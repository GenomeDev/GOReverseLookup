from collections import defaultdict
from typing import Union

from ..util.CacheUtil import Cacher

import requests
from requests.adapters import HTTPAdapter, Retry

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class gProfiler:
    def __init__(self) -> None:
        # Set up a retrying session
        retry_strategy = Retry(
            total=3, status_forcelist=[429, 500, 502, 503, 504], backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session

    # implementation as per https://biit.cs.ut.ee/gprofiler/convert
    def convert_ids(
        self,
        source_ids: Union[str, list[str], set[str]],
        taxon: str,
        target_namespace: str = "ensg",
    ) -> dict[str, list[str]]:
        """_summary_

        Args:
            source_ids (Union[str, list[str], set[str]]): _description_
            taxon (str): _description_
            target_namespace (str, optional): _description_. Defaults to "ensg".

        Raises:
            TypeError: _description_

        Returns:
            dict[str, list[str]]: _description_
        """
        USE_UNIPROT_IDMAP_NOTATION = True
        # example notation:
        # {
        #   'results': [{'from': 'P15692', 'to': 'ENSG00000112715.26'}, {'from': 'P16612', 'to': 'ENSRNOG00000019598'}, {'from': 'P40763', 'to': 'ENSG00000168610.16'}, {'from': 'P42227', 'to': 'ENSMUSG00000004040'}, {'from': 'P52631', 'to': 'ENSRNOG00000019742'}],
        #   'failedIds': ['O73682']
        # }

        if not isinstance(taxon, str):
            raise TypeError("taxons must be str")
        if isinstance(source_ids, str):
            source_ids_list = [source_ids]
        if isinstance(source_ids, set):
            source_ids_list = list(source_ids)
        if isinstance(source_ids, list):
            source_ids_list = source_ids

        converted_taxon = gProfilerUtil.NCBITaxon_to_gProfiler(taxon)
        if not converted_taxon:
            logger.warning(f"Failed to convert taxon {taxon}!")
            return {}
        namespace = target_namespace

        r = self.s.post(
            url="https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
            json={
                "organism": converted_taxon,
                "target": namespace,
                "query": source_ids_list,
            },
        )

        converted_ids = defaultdict(
            list, {k: [] for k in source_ids_list}
        )  # initialise with keys
        result: list[dict] = r.json()["result"]

        failedIds = []
        for entry in result:
            entry_source_id = entry["incoming"]
            if entry["converted"] not in ["N/A", "None", None]:
                converted_ids[entry_source_id].append(entry["converted"])
            else:
                failedIds.append(entry_source_id)

        if USE_UNIPROT_IDMAP_NOTATION == False:
            return converted_ids

        # return idmap as uniprot notation
        idmap_uniprot_notation_results = []
        for from_id, to_id in converted_ids.items():
            res = {"from": from_id, "to": to_id}
            idmap_uniprot_notation_results.append(res)

        idmap_uniprot_notation = {
            "results": idmap_uniprot_notation_results,
            "failedIds": failedIds,
        }
        return idmap_uniprot_notation

    def find_orthologs(
        self,
        source_ids: Union[str, list[str], set[str]],
        source_taxon: str,
        target_taxon: str = "9606",
    ) -> dict[str, list[str]]:
        """
        Maps the entire source_ids to respective 'target_taxon' Ensembl (ENS) gene ids.

        Args:
            source_ids (Union[str, list[str], set[str]]): Source ids which should be mapped to ENS ids. Example list is: ["ZDB-GENE-100922-278", "ZDB-GENE-131121-482", ...]
            source_taxon (str): The NCBITaxon number of the taxon corresponding to source ids. For example, if source ids are for Danio Rerio (Zebrafish), then source taxon should be "7955".
            target_taxon (str, optional): Defaults to "9606". Determines the final ENS id mapping (for example, if target taxon is 9606 (Homo Sapiens),
                                          then the mappings will be in the ENSG format. If the target taxon is Danio Rario, then the id mappings will be in the ENSDARG format)

        Returns:
            A dictionary mapping source_ids to respective ENS ids for target taxon.

            Example:
            {
                'ZDB-GENE-100922-278': ['ENSG00000160951'],
                'ZDB-GENE-131121-482': [],
                ...
            }
        """
        logger.info(f"[gProfiler::find_orthologs]: len source_ids = {len(source_ids)}, source_taxon = {source_taxon}, target_taxon = {target_taxon}")
        if source_ids == [] or source_ids == "":
            return None

        if ":" in source_taxon:
            source_taxon = source_taxon.split(":")[1]

        if not isinstance(source_taxon, str) and not isinstance(target_taxon, str):
            raise TypeError("taxons must be str")

        if isinstance(source_ids, str):
            source_ids_list = [source_ids]
        if isinstance(source_ids, set):
            source_ids_list = list(source_ids)
        if isinstance(source_ids, list):
            source_ids_list = source_ids

        target_ids = defaultdict(
            list, {k: [] for k in source_ids_list}
        )  # initialise with keys

        source_taxon = gProfilerUtil.NCBITaxon_to_gProfiler(source_taxon)
        target_taxon = gProfilerUtil.NCBITaxon_to_gProfiler(target_taxon)
        if not source_taxon or not target_taxon:
            return {}

        # cache previous data -> query only NEW ids
        new_input_ids = []  # input ids that were not yet queried
        cached_results = {}
        for source_id in source_ids:
            gOrth_data_key = f"[{self.__class__.__name__}][{self.find_orthologs.__name__}][id={source_id},source_taxon={source_taxon},target_taxon={target_taxon}]"
            previous_result = Cacher.get_data("gprofiler", gOrth_data_key)
            if previous_result is not None:
                cached_results[source_id] = previous_result
            else:
                new_input_ids.append(source_id)

        logger.info(f"Number of ids to be queried: {len(new_input_ids)}, source_taxon = {source_taxon}, target_taxon = {target_taxon}")
        if new_input_ids != []:
            """ # old link - problems?
            r = self.s.post(
                url="https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/orth/orth/",
                json={
                    "organism": source_taxon,
                    "target": target_taxon,
                    "query": new_input_ids,
                },
            )
            """ 
            r = self.s.post(
                url="https://biit.cs.ut.ee/gprofiler/api/orth/orth/",
                json={
                    "organism": source_taxon,
                    "target": target_taxon,
                    "query": new_input_ids,
                },
            )
            r.raise_for_status()
            result: list[dict] = r.json()["result"]

            # parse web-queried ids
            num_no_orthologs = 0
            num_orthologs = 0
            for entry in result:
                entry_source_id = entry["incoming"]
                gOrth_data_key = f"[{self.__class__.__name__}][{self.find_orthologs.__name__}][id={entry_source_id},source_taxon={source_taxon},target_taxon={target_taxon}]"
                if entry["ortholog_ensg"] not in ["N/A", "None", None]:
                    target_ids[entry_source_id].append(entry["ortholog_ensg"])
                    num_orthologs += 1
                    # don't cache here, since multiple ensgs can be appended!
                else:
                    # store no orthologs here!
                    num_no_orthologs += 1
                    Cacher.store_data("gprofiler", gOrth_data_key, "none")

            # eliminate double ids
            for entry_source_id, orthologs in target_ids.items():
                orthologs = set(orthologs)
                orthologs = list(orthologs)
                target_ids[entry_source_id] = orthologs

            # store web-queried ids in cache
            for entry_source_id, orthologs in target_ids.items():
                gOrth_data_key = f"[{self.__class__.__name__}][{self.find_orthologs.__name__}][id={entry_source_id},source_taxon={source_taxon},target_taxon={target_taxon}]"
                if orthologs == []:
                    Cacher.store_data("gprofiler", gOrth_data_key, "none")
                else:
                    Cacher.store_data("gprofiler", gOrth_data_key, orthologs)

            logger.debug(f"gProfiler orth query ({source_taxon} -> {target_taxon}): {len(new_input_ids)} input ids -> {num_orthologs} found, {num_no_orthologs} not found.")

        # parse cached ids
        for source_id, cached_value in cached_results.items():
            if cached_value != "none":
                if isinstance(cached_value, list):
                    # target_ids[source_id].append(*cached_value)
                    if len(cached_value) > 1:  # delete duplicates
                        cached_value = set(cached_value)
                        cached_value = list(cached_value)
                    target_ids[source_id] = [
                        *target_ids[source_id],
                        *cached_value,
                    ]  # join both lists
                else:
                    target_ids[source_id].append(cached_value)

        return target_ids


class gProfilerUtil:
    def __init__():
        pass

    @classmethod
    def NCBITaxon_to_gProfiler(cls, taxon: Union[str, int]):
        """
        Converts an NCBI-type taxon to a respective GProfiler taxon.
        Note: gprofiler - https://biit.cs.ut.ee/gprofiler/

        Args:
            taxon (str or int): a full NCBIType-taxon (NCBITaxon:xxxx) or an integer representing the NCBITaxon id (xxxx)

        Returns:
            _type_: _description_
        """
        if isinstance(taxon, str):
            if ":" in taxon:
                taxon = taxon.split(":")[1]

        url = "https://biit.cs.ut.ee/gprofiler/api/util/organisms_list"
        prev_response = Cacher.get_data("gprofiler", url)
        if prev_response is None:
            r = requests.get(url)
            results = r.json()
            Cacher.store_data("gprofiler", data_key=url, data_value=results)
        else:
            results = prev_response
            
        taxon_equivalents = {}
        for r in results:
            taxon_equivalents[r["taxonomy_id"]] = r["id"]
        return taxon_equivalents.get(str(taxon), None)

from collections import defaultdict
from typing import Union

import requests
from requests.adapters import HTTPAdapter, Retry


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
        
        if not isinstance(taxon, str):
            raise TypeError("taxons must be str")
        if isinstance(source_ids, str):
            source_ids_list = [source_ids]
        if isinstance(source_ids, set):
            source_ids_list = list(source_ids)
        if isinstance(source_ids, list):
            source_ids_list = source_ids

        converted_taxon = NCBITaxon_to_gProfiler(taxon)
        if not converted_taxon:
            return {}
        namespace = target_namespace

        r = self.s.post(
            url="https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
            json={
                "organism": taxon,
                "target": namespace,
                "query": source_ids_list,
            },
        )

        converted_ids = defaultdict(
            list, {k: [] for k in source_ids_list}
        )  # initialise with keys
        result: list[dict] = r.json()["result"]

        for entry in result:
            entry_source_id = entry["incoming"]
            if entry["converted"] not in ["N/A", "None", None]:
                converted_ids[entry_source_id].append(entry["converted"])

        return converted_ids


def NCBITaxon_to_gProfiler(taxon):
    """_summary_

    Args:
        taxon (_type_): _description_

    Returns:
        _type_: _description_
    """
    r = requests.get("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list")
    taxon_equivalents = {}
    results = r.json()
    for r in results:
        taxon_equivalents[r["taxonomy_id"]] = r["id"]
    return taxon_equivalents.get(str(taxon), None)

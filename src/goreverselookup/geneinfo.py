from collections import defaultdict
from typing import Union

import requests

from .utils import NCBITaxon_to_gProfiler


def gConvert(ids: list[str], taxon, namespace: str) -> dict[str, list[str]]:
    """_summary_

    Args:
        ids (list[str]): _description_
        taxon (_type_): _description_
        namespace (str): _description_

    Returns:
        dict[str, list[str]]: _description_
    """
    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
        json={
            "organism": taxon,
            "target": namespace,
            "query": ids,
        },
    )

    converted_ids = defaultdict(list, {k: [] for k in ids})  # initialise with keys
    result: list[dict] = r.json()["result"]

    for entry in result:
        entry_source_id = entry["incoming"]
        if entry["converted"] not in ["N/A", "None", None]:
            converted_ids[entry_source_id].append(entry["converted"])
    return converted_ids


def convert_ids(
    source_ids: Union[str, list[str], set[str]],
    taxon: str,
    target_namespace: str = "ensg",
    database: str = "gConvert",
):
    """_summary_

    Args:
        source_ids (Union[str, list[str], set[str]]): _description_
        taxon (str): _description_
        target_namespace (str, optional): _description_. Defaults to "ensg".
        database (str, optional): _description_. Defaults to "gConvert".

    Raises:
        TypeError: _description_

    Returns:
        _type_: _description_
    """
    if not isinstance(taxon, str):
        raise TypeError("taxons must be str")
    if isinstance(source_ids, str):
        source_ids_list = [source_ids]
    if isinstance(source_ids, set):
        source_ids_list = list(source_ids)
    if isinstance(source_ids, list):
        source_ids_list = source_ids

    if database == "gConvert":
        converted_taxon = NCBITaxon_to_gProfiler(taxon)
        if not converted_taxon:
            return {}
        namespace = target_namespace
        converted_ids = gConvert(source_ids_list, converted_taxon, namespace)
    else:
        return ValueError(f"database {database} is not available.")

    return converted_ids

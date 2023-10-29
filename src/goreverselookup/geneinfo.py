from collections import defaultdict
from typing import Union

import requests
import json
import aiohttp
import asyncio

from .utils import NCBITaxon_to_gProfiler, NCBITaxon_to_ensembl


def gConvert(ids: list[str], taxon, namespace: str) -> dict[str, set[str]]:
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

    converted_ids: dict[str, set[str]] = defaultdict(
        set, {k: set() for k in ids}
    )  # initialise with keys
    result: list[dict] = r.json()["result"]

    for entry in result:
        entry_source_id = entry["incoming"]
        if entry["converted"] not in ["N/A", "None", None]:
            converted_ids[entry_source_id].add(entry["converted"])
    return converted_ids


async def async_ensembl_xrefs(ids: list[str], taxon) -> dict[str, set[str]]:
    """_summary_

    Args:
        ids (list[str]): _description_
        taxon (_type_): _description_

    Returns:
        dict[str, set[str]]: _description_
    """
    list_of_urls = []
    for id in ids:
        symbol = id.split(":", 1)[1] if len(id.split(":", 1)) > 1 else id
        url = f"https://rest.ensembl.org/xrefs/symbol/{taxon}/{symbol}?content-type=application/json"
        list_of_urls.append(url)

    async def _semaphore_get_request(url, session, sem):
        async with sem:
            async with session.get(url) as response:
                data = await response.json()
            return data

    async with aiohttp.ClientSession() as session:
        semaphore = asyncio.Semaphore(100)  # max tries per second

        # Create and start tasks with rate limiting
        tasks = [
            _semaphore_get_request(url, session, semaphore) for url in list_of_urls
        ]

        results = await asyncio.gather(*tasks)

    converted_ids: dict[str, set[str]] = defaultdict(
        set, {k: set() for k in ids}
    )  # initialise with keys

    for id, entry in zip(ids, results):
        if not entry:
            pass
        for converted_dict in entry:
            if not converted_dict["type"] == "gene":
                pass
            converted_ids[id].add(converted_dict["id"])
    return converted_ids


def ensembl_xrefs(ids: list[str], taxon) -> dict[str, set[str]]:
    return asyncio.run(async_ensembl_xrefs(ids, taxon))

def convert_ids(
    source_ids: Union[str, list[str], set[str]],
    taxon: str,
    target_namespace: str = "ensg",
    database: str = "gConvert",
) -> dict[str, set[str]]:
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
        namespace = target_namespace  # in future maybe filter?
        converted_ids = gConvert(source_ids_list, converted_taxon, namespace)
    elif database == "ensembl":
        converted_taxon = NCBITaxon_to_ensembl(taxon)
        if not converted_taxon:
            return {}
        if target_namespace != "ensg":
            raise ValueError("database ensembl can only provide ensg ids")
        converted_ids = ensembl_xrefs(source_ids_list, converted_taxon)
    else:
        raise ValueError(f"database {database} is not available.")

    return converted_ids

def ensembl_lookup(ids: list[str], taxon) -> dict[str, dict[str, str]]:
    """TODO:_summary_

    Args:
        ids (list[str]): _description_
        taxon (_type_): _description_

    Returns:
        dict[str, dict[str, str]]: _description_
    """

    # Split the IDs into chunks of 1000 or less (ensembl api limit)
    chunk_size = 1000
    ids_chunks = [ids[i : i + chunk_size] for i in range(0, len(ids), chunk_size)]

    results = {}
    for chunk in ids_chunks:
        r = requests.post(
            url=f"https://rest.ensembl.org/lookup/symbol/{taxon}",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            data=json.dumps({"symbols": chunk}),
        )
        if r:
            results.update(r.json())

    ids_dict: dict[str, dict[str, str]] = defaultdict(dict, {k: {} for k in ids})
    for id, result in results.items():
        ids_dict[id]["ensg"] = result["id"]
        ids_dict[id]["genename"] = result["display_name"]
        ids_dict[id]["description"] = result["description"]

    return ids_dict

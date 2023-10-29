from collections import defaultdict
from typing import Union

import requests
import aiohttp
import asyncio

from .utils import NCBITaxon_to_gProfiler, NCBITaxon_to_ensembl


def gOrth(
    source_ids: list[str], source_taxon: str, target_taxon: str
) -> dict[str, set[str]]:
    """_summary_

    Args:
        source_ids (list): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_
    """
    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/orth/orth/",
        json={
            "organism": source_taxon,
            "target": target_taxon,
            "query": source_ids,
        },
    )

    target_ids = defaultdict(
        set, {k: set() for k in source_ids}
    )  # initialise with keys
    result: list[dict] = r.json()["result"]
    for entry in result:
        entry_source_id = entry["incoming"]
        if entry["ortholog_ensg"] not in ["N/A", "None", None]:
            target_ids[entry_source_id].add(entry["ortholog_ensg"])
    return target_ids


async def async_ensembl_orthologs(
    source_ids: list[str], source_taxon: str, target_taxon: str
) -> dict[str, set[str]]:
    """_summary_

    Args:
        source_ids (list[str]): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_

    Returns:
        dict[str, list[str]]: _description_
    """

    list_of_urls = []
    for source_id in source_ids:
        symbol = (
            source_id.split(":", 1)[1]
            if len(source_id.split(":", 1)) > 1
            else source_id
        )
        url = f"https://rest.ensembl.org/homology/symbol/{source_taxon}/{symbol}?target_species={target_taxon};type=orthologues;sequence=none;content-type=application/json"
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

    target_ids = defaultdict(
        set, {k: set() for k in source_ids}
    )  # initialise with keys
    for source_id, result in zip(source_ids, results):
        if result == [] or "error" in result:
            pass

        homologies_list_of_dicts = result.get("data", [{}])[0].get(
            "homologies"
        )  # is result["data"][0]["homologies"]
        if not homologies_list_of_dicts:
            pass

        best_ortholog = max(
            homologies_list_of_dicts, key=lambda a: a["target"]["perc_id"]
        )

        target_ids[source_id].add(best_ortholog["target"]["id"])

    return target_ids


def ensembl_orthologs(
    source_ids: list[str], source_taxon: str, target_taxon: str
) -> dict[str, set[str]]:
    return asyncio.run(async_ensembl_orthologs(source_ids, source_taxon, target_taxon))


def find_orthologs(
    source_ids: Union[str, list[str], set[str]],
    source_taxon: str,
    target_taxon: str = "9606",
    database: str = "gOrth",
) -> dict[str, set[str]]:
    """_summary_

    Args:
        source_ids (Union[str, list[str], set[str]]): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_
        database (str): _description_

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: _description_
    """
    if not isinstance(source_taxon, str) and not isinstance(target_taxon, str):
        raise TypeError("taxons must be str")

    if isinstance(source_ids, str):
        source_ids_list = [source_ids]
    if isinstance(source_ids, set):
        source_ids_list = list(source_ids)
    if isinstance(source_ids, list):
        source_ids_list = source_ids

    if database == "gOrth":
        source_taxon = NCBITaxon_to_gProfiler(source_taxon)
        target_taxon = NCBITaxon_to_gProfiler(target_taxon)
        if not source_taxon or not target_taxon:
            return {}
        target_ids_dict = gOrth(source_ids_list, source_taxon, target_taxon)
    elif database == "ensembl":
        source_taxon = NCBITaxon_to_ensembl(source_taxon)
        target_taxon = NCBITaxon_to_ensembl(target_taxon)
        if not source_taxon or not target_taxon:
            return {}
        target_ids_dict = ensembl_orthologs(source_ids_list, source_taxon, target_taxon)
    elif database == "local_files":
        raise NotImplementedError
    else:
        ValueError(
            f"database {database} is not available as a source of ortholog information"
        )

    return target_ids_dict

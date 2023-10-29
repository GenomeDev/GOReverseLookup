import requests


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


def NCBITaxon_to_ensembl(taxon):
    """_summary_

    Args:
        taxon (_type_): _description_
    """
    r = requests.get(
        "https://rest.ensembl.org/info/species?content-type=application/json"
    )
    taxon_equivalents = {}
    results = r.json()
    for r in results["species"]:
        taxon_equivalents[r["taxon_id"]] = r["name"]
    return taxon_equivalents.get(str(taxon), None)

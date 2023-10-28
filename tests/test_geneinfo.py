from goreverselookup.geneinfo import gConvert


def test_gConvert():
    results = gConvert(["ZDB-GENE-021119-1"], "drerio", "ensg")
    assert len(results["ZDB-GENE-021119-1"]) == 2

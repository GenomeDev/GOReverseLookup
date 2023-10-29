import pytest
from goreverselookup.geneinfo import gConvert, async_ensembl_xrefs


def test_gConvert():
    results = gConvert(["ZDB-GENE-021119-1"], "drerio", "ensg")
    assert results["ZDB-GENE-021119-1"] == {"ENSDARG00000032131", "ENSDARG00000110679"}


@pytest.mark.asyncio
async def test_async_ensembl_xrefs():
    results = await async_ensembl_xrefs(["ZDB-GENE-021119-1"], "danio_rerio")
    assert results["ZDB-GENE-021119-1"] == {"ENSDARG00000032131", "ENSDARG00000110679"}

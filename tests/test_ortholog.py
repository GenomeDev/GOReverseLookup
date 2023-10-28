import pytest

from goreverselookup.ortholog import find_orthologs, gOrth, ensemblOrthologs


def test_gOrth():
    result = gOrth(["ZDB-GENE-040912-6", "ZDB-GENE-170217-1"], "drerio", "hsapiens")
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000005421", "ENSG00000105852", "ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
    }


def test_missing_gOrth():
    result = gOrth(["A0A087WV62"], "drerio", "hsapiens")
    assert result == {
        "A0A087WV62": [],
    }


def test_find_orthologs():
    result = find_orthologs(
        ["ZDB-GENE-040912-6", "ZDB-GENE-170217-1", "ZDB-GENE-021119-1"],
        "7955",
        "9606",
        database="gOrth",
    )
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000005421", "ENSG00000105852", "ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
        "ZDB-GENE-021119-1": ["ENSG00000151577"],
    }


@pytest.mark.parametrize("db", ["gOrth", "ensembl"])
def test_multiple_same_incoming_find_orthologs(db):
    result = find_orthologs(
        ["ZDB-GENE-170217-1", "ZDB-GENE-040912-6"], "7955", "9606", database=db
    )
    assert "ENSG00000105854" in result["ZDB-GENE-040912-6"]
    assert "ENSG00000168938" in result["ZDB-GENE-170217-1"]

@pytest.mark.asyncio
async def test_ensemblOrtholog():
    result = await ensemblOrthologs(["ZDB-GENE-040912-6", "ZDB-GENE-170217-1"], "danio_rerio", "homo_sapiens")
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
    }
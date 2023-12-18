import pytest

from goreverselookup.web_apis.gProfilerApi import gProfiler
from goreverselookup import gProfilerUtil


@pytest.fixture
def gProfiler_class():
    return gProfiler()


def test_convert_ids(gProfiler_class):
    results = gProfiler_class.convert_ids(["ZDB-GENE-021119-1"], "7955", "ensg")
    assert len(results["ZDB-GENE-021119-1"]) == 2


def test_NCBITaxon_to_gProfiler():
    assert gProfilerUtil.NCBITaxon_to_gProfiler("7955") == "drerio"


def test_find_orthologs(gProfiler_class):
    result = gProfiler_class.find_orthologs(
        ["ZDB-GENE-040912-6", "ZDB-GENE-170217-1", "ZDB-GENE-021119-1"],
        "7955",
        "9606",
    )
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000005421", "ENSG00000105852", "ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
        "ZDB-GENE-021119-1": ["ENSG00000151577"],
    }

import pytest

from goreverselookup.web_apis.gProfilerApi import NCBITaxon_to_gProfiler, gProfiler


@pytest.fixture
def gProfiler_class():
    return gProfiler()


def test_convert_ids(gProfiler_class):
    results = gProfiler_class.convert_ids(["ZDB-GENE-021119-1"], "7955", "ensg")
    assert len(results["ZDB-GENE-021119-1"]) == 2


def test_NCBITaxon_to_gProfiler():
    assert NCBITaxon_to_gProfiler("7955") == "drerio"

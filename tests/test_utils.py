import pytest

from goreverselookup.utils import NCBITaxon_to_ensembl, NCBITaxon_to_gProfiler


def test_NCBITaxon_to_ensembl():
    assert NCBITaxon_to_ensembl("9606") == "homo_sapiens"


def test_NCBITaxon_to_gProfiler():
    assert NCBITaxon_to_gProfiler("9606") == "hsapiens"

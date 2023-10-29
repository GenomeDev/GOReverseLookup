"""
Here should be all fixtures, so that they can be reused by tests.
A fixture is a function which prepares fresh data (object,list,string...) for a test.
"""

import pytest

from goreverselookup.web_apis.EnsemblApi import EnsemblApi
from goreverselookup import Cacher
import aiohttp


@pytest.fixture
def cacher(tmp_path):
    # this creates unique cache in a temp folder
    return Cacher.init(cache_dir=tmp_path)


@pytest.fixture
def ensemblapi(cacher) -> EnsemblApi:
    # this creates unique ensembl api object with unique cache for each test
    cacher
    return EnsemblApi()


@pytest.fixture
async def aiohttp_session():
    connector = aiohttp.TCPConnector(limit=100, limit_per_host=100)
    async with aiohttp.ClientSession(connector=connector) as session:
        yield session


ortholog_tuple_list = [
    ("MGI:88190", "ENSG00000157764"),
    ("ZFIN:ZDB-GENE-040805-1", "ENSG00000157764"),
    ("RGD:619908", "ENSG00000157764"),
    ("Xenbase:XB-GENE-1014202", "ENSG00000157764"),
]


@pytest.fixture(params=ortholog_tuple_list)
def ortholog_tuple_parametrized(request):
    # it is parametrized -> a test which uses it will run
    # 4 times each time with different parameters
    return request.param[0], request.param[1]


ensembl_data_list = [
    (
        "ENSG00000157764",
        {
            "ensg_id": "ENSG00000157764",
            "genename": "BRAF",
            "description": "B-Raf proto-oncogene, serine/threonine kinase",
            "enst_id": "ENST00000644969",
            "refseq_nt_id": "NM_001374258.1",
            "uniprot_id": "",
        },
    ),
    (
        "MGI:88190",
        {
            "ensg_id": "ENSMUSG00000002413",
            "genename": "Braf",
            "description": "Braf transforming gene",
            "enst_id": "ENSMUST00000002487",
            "refseq_nt_id": None,
            "uniprot_id": "P28028",
        },
    ),
]


@pytest.fixture(params=ensembl_data_list)
def ensembl_data_parametrized(request):
    return request.param[0], request.param[1]

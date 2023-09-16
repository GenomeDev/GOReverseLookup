import pytest

from goreverselookup.web_apis.EnsemblApi import EnsemblApi
from goreverselookup import Cacher

@pytest.fixture
def cacher(tmp_path):
    #this creates unique cache in a temp folder
    return Cacher.init(cache_dir=tmp_path)

@pytest.fixture
def ensemblapi(cacher):
    #this creates unique ensembl api object with unique cache for each test
    cacher
    return EnsemblApi()

@pytest.mark.usefixtures("ensemblapi")
class TestEnsemblAPI():
    @pytest.mark.parametrize("source_id, expected", [("MGI:88190", "ENSG00000157764")])
    def test_get_human_ortholog(self, ensemblapi, source_id, expected):
        ortholog_id = ensemblapi.get_human_ortholog(source_id)
        assert ortholog_id == expected
        
        
    

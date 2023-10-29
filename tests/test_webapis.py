import pytest
import re
import logging

from goreverselookup.web_apis.EnsemblApi import EnsemblApi


class TestEnsemblAPI:
    def test_get_human_ortholog(
        self, caplog, ensemblapi: EnsemblApi, ortholog_tuple_parametrized
    ):
        # test for each of the supported species
        source_id, expected = ortholog_tuple_parametrized
        ortholog_id = ensemblapi.get_human_ortholog(source_id)
        assert ortholog_id == expected

    def test_get_human_ortholog_cache(self, caplog, ensemblapi: EnsemblApi):
        # test if cache works (by listening to logger output)
        caplog.set_level(logging.DEBUG)
        ensemblapi.get_human_ortholog("MGI:88190")
        ensemblapi.get_human_ortholog("MGI:88190")  # run it twice
        assert "cached" in caplog.records[-1].message

    @pytest.mark.skip
    @pytest.mark.asyncio
    async def test_get_human_ortholog_async(
        self, ensemblapi: EnsemblApi, aiohttp_session, ortholog_tuple
    ):
        # test for each of the supported species
        for pair in ortholog_tuple:
            source_id, expected = pair
            ortholog_id = await ensemblapi.get_human_ortholog_async(
                source_id, session=aiohttp_session
            )
            assert ortholog_id == expected

    def test_get_sequence(self, ensemblapi: EnsemblApi):
        ensg_id = "ENST00000646891"
        no_type_arg_seq = ensemblapi.get_sequence(ensg_id)
        assert re.match(
            "^[ATGC]*$", no_type_arg_seq
        )  # assert that the string only has A,T,C and G in it.
        cdna_seq = ensemblapi.get_sequence(ensg_id, sequence_type="cdna")
        assert cdna_seq == no_type_arg_seq  # assert that default "cdna" works
        cds_seq = ensemblapi.get_sequence(ensg_id, sequence_type="cds")
        assert re.match(
            "^[ATGC]*$", cds_seq
        )  # assert that the string only has A,T,C and G in it.
        assert cds_seq != cdna_seq  # assert that cds is not the same as cdna

    def test_get_info_bugfix(self, ensemblapi: EnsemblApi):
        # test the bugfix for "error" in id
        assert (
            ensemblapi.get_info(
                "RgdError_No-human-ortholog-found:product_id=RGD:1359312"
            )
            == {}
        )

    def test_get_info(self, ensemblapi: EnsemblApi, ensembl_data_parametrized):
        source_id, expected = ensembl_data_parametrized
        info = ensemblapi.get_info(source_id)
        assert info == expected

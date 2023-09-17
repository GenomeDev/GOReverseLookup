from typing import Optional, List
import asyncio
import aiohttp

from .ModelSettings import ModelSettings
from ..parse.OrthologParsers import HumanOrthologFinder
from ..parse.GOAFParser import GOAnnotationsFile
from ..web_apis.EnsemblApi import EnsemblApi
from ..web_apis.UniprotApi import UniProtApi

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class Product:
    def __init__(
        self,
        id_synonyms: List[str],
        genename: str = None,
        uniprot_id: str = None,
        description: str = None,
        ensg_id: str = None,
        enst_id: str = None,
        refseq_nt_id: str = None,
        mRNA: str = None,
        scores: dict = None,
        had_orthologs_computed: bool = False,
        had_fetch_info_computed: bool = False,
    ):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id_synonyms (str): The list of ID of the product and synonyms. -> after ortholog translation it turns out that some products are the same. Example: RGD:3774, Xenbase:XB-GENE-5818802, UniProtKB:Q9NTG7
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            ensg_id (str): Ensembl gene ID (MAIN).
            enst_id (str): Ensembl transcript ID.
            refseq_nt_id (str): Refseq (reference sequence) transcript ID.
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
            had_orthologs_computed (bool): If this Product instance has had the fetch_ortholog function called already.
            had_fetch_info_computed (bool): If this Product instance has had the fetch_info function called already.
        """
        self.id_synonyms = id_synonyms
        self.genename = genename  # NOTE: genename indicates a successful ortholog fetch operation !!!
        self.description = description
        self.uniprot_id = uniprot_id
        self.ensg_id = ensg_id
        self.enst_id = enst_id
        self.refseq_nt_id = refseq_nt_id
        self.mRNA = mRNA
        self.scores = {} if scores is None else scores.copy()
        self.had_orthologs_computed = had_orthologs_computed
        self.had_fetch_info_computed = had_fetch_info_computed
        self._d_offline_online_ortholog_mismatch = False  # if fetch_ortholog is queried with _d_compare_goaf set to True, this variable will be set to True if there is a mismatch in the gene names returned from the online and offline query algorithms.
        self._d_offline_online_ortholog_mismatch_values = ""

        # see if UniProtKB id is already in id_synonyms:
        for id_syn in self.id_synonyms:
            if "UniProt" in id_syn:
                self.uniprot_id = id_syn

    def fetch_ortholog(
        self,
        human_ortholog_finder: Optional[HumanOrthologFinder] = None,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        goaf: Optional[GOAnnotationsFile] = None,
        prefer_goaf=False,
        _d_compare_goaf=False,
        model_settings: Optional[ModelSettings] = None,
    ) -> None:
        """
        Fetches the ortholog for this product. If the ortholog query was successful, then self.genename is updated to the correct human ortholog gene name.
        Additionally, during the course of fetch_ortholog, ensembl_api.get_info may be called - if this happens, then the values description, ensg_id, enst_id, refseq_nt_id, uniprot_id are
        also filled out for this Product from the ensembl_api.get_info return value.

        Parameters:
          - (HumanOrthologFinder) human_ortholog_finder
          - (UniProtAPI) uniprot_api
          - (EnsemblAPI) ensembl_api
          - (GOAnnotationsFile) goaf
          - (bool) prefer_goaf: see explanation in the Algorithm section
          - (bool) _d_compare_goaf: if true, will attempt ortholog search both from offline and online algorithms and report if the results are the same
          - (ModelSettings) model_settings: isn't implemented, as this function is outdated. Use this function to check _d_compare_goaf!

        Algorithm:
            If there is only one id_synonym (eg. UniProtKB:Q9NTG7) and that id_synonym is of type UniProtKB, then
            UniProtAPI is used to obtained information about this gene. A successful query returns a dictionary, which also
            contains the genename field (which updates self.genename to the queried genename)

            If there is only one id_synonym and it is not of type UniProtKB, then HumanOrthologFinder is used to attempt a file-based
            search for the ortholog (files from all the third party databases are used).

            The user also has an option to supply a GO Annotations File and direct the program to first browse the GOAF and the 3rd party
            database files for orthologs ("offline" approach) using the prefer_goaf parameter. By default, if a GOAF is provided, it will be preferably used.

            If the file-based search doesn't work, then EnsemblAPI is used as a final attempt to find a human ortholog. The first call (ensembl_api.get_human_ortholog)
            returns an ensg_id, which is then used in another call to ensembl_api.get_info in order to obtain the gene name from the ensg_id.

            TODO: If there are multiple id_synonym(s), currently only the first is browsed. Implement logic for many id_synonyms / check if there are any products with multiple id synonyms.
        """
        DROP_MIRNA_FROM_ENSEMBL_QUERY = True  # returns None if Ensembl query returns a miRNA (MIRxxx) as the gene name.

        if not human_ortholog_finder:
            human_ortholog_finder = HumanOrthologFinder()
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        # *** offline (GOAF) and 3rd-party-database-file based analysis ***
        offline_queried_ortholog = None
        if prefer_goaf is True or _d_compare_goaf is True:
            if len(self.id_synonyms) == 1 and "UniProtKB" in self.id_synonyms[0]:
                # find the DB Object Symbol in the GOAF. This is the third line element. Example: UniProtKB	Q8NI77	KIF18A	located_in	GO:0005737	PMID:18680169	IDA		C	Kinesin-like protein KIF18A	KIF18A|OK/SW-cl.108	protein	taxon:9606	20090818	UniProt -> KIF18A
                if goaf is not None:
                    self.genename = goaf.get_uniprotkb_genename(self.id_synonyms[0])
                else:
                    logger.warning(
                        "GOAF wasn't supplied as parameter to the"
                        " (Product).fetch_ortholog function!"
                    )
            elif len(self.id_synonyms) == 1 and "UniProtKB" not in self.id_synonyms[0]:
                # do a file-based ortholog search using HumanOrthologFinder
                human_ortholog_gene_id = human_ortholog_finder.find_human_ortholog(
                    self.id_synonyms[0]
                )
                offline_queried_ortholog = human_ortholog_gene_id  # this is used for acceleration so as not to repeat find_human_ortholog in the online algorithm section
                if human_ortholog_gene_id is not None:
                    self.genename = human_ortholog_gene_id

        # *** online and 3rd-party-database-file based analysis ***
        if _d_compare_goaf is True or prefer_goaf is False:
            if len(self.id_synonyms) == 1 and "UniProtKB" in self.id_synonyms[0]:
                if self.uniprot_id is None:
                    # 14.08.2023: replaced online uniprot info query with goaf.get_uniprotkb_genename, as it is more successful and does the same as the uniprot query
                    # online uniprot info query is performed only for debugging purposes with _d_compare_goaf
                    if _d_compare_goaf is True:
                        info_dict = uniprot_api.get_uniprot_info(
                            self.id_synonyms[0]
                        )  # bugfix
                    else:
                        info_dict = {
                            "genename": goaf.get_uniprotkb_genename(self.id_synonyms[0])
                        }
                else:  # self.uniprot_id exists
                    if _d_compare_goaf is True:
                        info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
                    else:
                        info_dict = {
                            "genename": goaf.get_uniprotkb_genename(self.uniprot_id)
                        }
                # if compare is set to True, then only log the comparison between
                if _d_compare_goaf is True:
                    if self.genename != info_dict.get("genename"):
                        logger.warning(
                            f"GOAF-obtained genename ({self.genename}) is not the same"
                            " as UniProtKB-obtained genename"
                            f" ({info_dict.get('genename')}) for {self.id_synonyms}"
                        )
                        self._d_offline_online_ortholog_mismatch = True
                        self._d_offline_online_ortholog_mismatch_values = (
                            f"[{self.id_synonyms[0]}]: online ="
                            f" {info_dict.get('genename')}, offline = {self.genename};"
                            " type = uniprot query"
                        )
                else:
                    self.genename = info_dict.get("genename")

            elif len(self.id_synonyms) == 1 and "UniProtKB" not in self.id_synonyms[0]:
                if (
                    offline_queried_ortholog is None
                ):  # if algorithm enters this section due to _d_compare_goaf == True, then this accelerates code, as it prevents double calculations
                    human_ortholog_gene_id = human_ortholog_finder.find_human_ortholog(
                        self.id_synonyms[0]
                    )  # file-based search; alternative spot for GOAF analysis
                else:
                    human_ortholog_gene_id = offline_queried_ortholog
                if (
                    human_ortholog_gene_id is None
                ):  # if file-based search finds no ortholog
                    logger.warning(
                        "human ortholog finder did not find ortholog for"
                        f" {self.id_synonyms[0]}"
                    )
                    human_ortholog_gene_ensg_id = ensembl_api.get_human_ortholog(
                        self.id_synonyms[0]
                    )  # attempt ensembl search
                    if human_ortholog_gene_ensg_id is not None:
                        enst_dict = ensembl_api.get_info(human_ortholog_gene_ensg_id)
                        human_ortholog_gene_id = enst_dict.get("genename")
                        if human_ortholog_gene_id is not None:
                            if (
                                DROP_MIRNA_FROM_ENSEMBL_QUERY is True
                                and "MIR" in human_ortholog_gene_id
                            ):
                                human_ortholog_gene_id = (
                                    None  # Ensembl query returned a miRNA, return None
                                )
                        if _d_compare_goaf is True:
                            if self.genename != human_ortholog_gene_id:
                                logger.warning(
                                    f"GOAF-obtained genename ({self.genename}) is not"
                                    " the same as Ensembl-obtained genename"
                                    f" ({human_ortholog_gene_id}) for"
                                    f" {self.id_synonyms}"
                                )
                                self._d_offline_online_ortholog_mismatch = True
                                self._d_offline_online_ortholog_mismatch_values = (
                                    f"[{self.id_synonyms[0]}]: online ="
                                    f" {human_ortholog_gene_id}, offline ="
                                    f" {self.genename}, type = ensembl query"
                                )
                        else:
                            self.genename = enst_dict.get("genename")
                            # update 19.08.2023: attempt to obtain as many values as possible for this Product already from
                            # the ortholog fetch to avoid duplicating requests with (EnsemblAPI).get_info
                            if self.ensg_id == "" or self.ensg_id is None:
                                self.ensg_id = enst_dict.get("ensg_id")
                            if self.description == "" or self.description is None:
                                self.description = enst_dict.get("description")
                            if self.enst_id == "" or self.enst_id is None:
                                self.enst_id = enst_dict.get("enst_id")
                            if self.refseq_nt_id == "" or self.refseq_nt_id is None:
                                self.refseq_nt_id == enst_dict.get("refseq_nt_id")
                            if self.uniprot_id == "" or self.uniprot_id is None:
                                uniprot_id = enst_dict.get("uniprot_id")
                                if uniprot_id is not None and uniprot_id != "":
                                    self.uniprot_id = enst_dict.get("uniprot_id")
                else:
                    if _d_compare_goaf is True:
                        if (
                            self.genename != human_ortholog_gene_id
                        ):  # with the current workflow, these will always be the same
                            logger.warning(
                                f"GOAF-obtained genename ({self.genename}) is not the"
                                " same as file-search-obtained-genename"
                                f" ({human_ortholog_gene_id}) for {self.id_synonyms}"
                            )
                    else:
                        self.genename = human_ortholog_gene_id

        self.had_orthologs_computed = True

    async def fetch_ortholog_async(
        self,
        session: aiohttp.ClientSession,
        goaf: GOAnnotationsFile,
        human_ortholog_finder: Optional[HumanOrthologFinder] = None,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        model_settings: Optional[ModelSettings] = None,
    ) -> None:
        """
        Fetches the ortholog for this product. If the ortholog query was successful, then self.genename is updated to the correct human ortholog gene name.
        Additionally, during the course of fetch_ortholog, ensembl_api.get_info may be called - if this happens, then the values description, ensg_id, enst_id, refseq_nt_id, uniprot_id are
        also filled out for this Product from the ensembl_api.get_info return value.

        Parameters:
          - (HumanOrthologFinder) human_ortholog_finder
          - (UniProtAPI) uniprot_api
          - (EnsemblAPI) ensembl_api
          - (GOAnnotationsFile) goaf
          - (bool) prefer_goaf: see explanation in the Algorithm section
          - (bool) _d_compare_goaf: if true, will attempt ortholog search both from offline and online algorithms and report if the results are the same
          - (ModelSettings) model_settings: the settings of the model. Currently, model_settings.uniprotkb_genename_online_query is used, which determines if gene name querying from a UniProtKB id is done via a web request or via GOAF

        Algorithm:
            If there is only one id_synonym (eg. UniProtKB:Q9NTG7) and that id_synonym is of type UniProtKB, then
            UniProtAPI is used to obtained information about this gene. A successful query returns a dictionary, which also
            contains the genename field (which updates self.genename to the queried genename)

            If there is only one id_synonym and it is not of type UniProtKB, then HumanOrthologFinder is used to attempt a file-based
            search for the ortholog (files from all the third party databases are used).

            The user also has an option to supply a GO Annotations File and direct the program to first browse the GOAF and the 3rd party
            database files for orthologs ("offline" approach) using the prefer_goaf parameter. By default, if a GOAF is provided, it will be preferably used.

            If the file-based search doesn't work, then EnsemblAPI is used as a final attempt to find a human ortholog. The first call (ensembl_api.get_human_ortholog)
            returns an ensg_id, which is then used in another call to ensembl_api.get_info in order to obtain the gene name from the ensg_id.

            TODO: If there are multiple id_synonym(s), currently only the first is browsed. Implement logic for many id_synonyms / check if there are any products with multiple id synonyms.
        """
        # logger.info(f"Async fetch orthologs for: {self.id_synonyms}")

        if not human_ortholog_finder:
            human_ortholog_finder = HumanOrthologFinder()
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        if len(self.id_synonyms) == 1 and "UniProtKB" in self.id_synonyms[0]:
            if self.uniprot_id is None or self.uniprot_id == "":
                # 14.08.2023: replaced online uniprot info query with goaf.get_uniprotkb_genename, as it is more successful and does the same as the uniprot query
                if (
                    model_settings is not None
                    and model_settings.uniprotkb_genename_online_query is True
                ):
                    info_dict = await uniprot_api.get_uniprot_info_async(
                        self.id_synonyms[0], session
                    )
                else:
                    info_dict = {
                        "genename": goaf.get_uniprotkb_genename(self.id_synonyms[0])
                    }
            else:
                if (
                    model_settings is not None
                    and model_settings.uniprotkb_genename_online_query is True
                ):
                    info_dict = await uniprot_api.get_uniprot_info_async(
                        self.uniprot_id, session
                    )
                else:
                    info_dict = {
                        "genename": goaf.get_uniprotkb_genename(self.uniprot_id)
                    }
            if info_dict is not None:
                self.genename = info_dict.get("genename")
        elif len(self.id_synonyms) == 1:
            human_ortholog_gene_id = (
                await human_ortholog_finder.find_human_ortholog_async(
                    self.id_synonyms[0]
                )
            )
            if human_ortholog_gene_id is None:
                logger.warning(
                    "human ortholog finder did not find ortholog for"
                    f" {self.id_synonyms[0]}"
                )
                human_ortholog_gene_ensg_id = (
                    await ensembl_api.get_human_ortholog_async(
                        self.id_synonyms[0], session
                    )
                )  # attempt ensembl search
                if human_ortholog_gene_ensg_id is not None:
                    enst_dict = await ensembl_api.get_info_async(
                        human_ortholog_gene_ensg_id, session
                    )
                    self.genename = enst_dict.get("genename")

                    # update 19.08.2023: attempt to obtain as many values as possible for this Product already from
                    # the ortholog fetch to avoid duplicating requests with (EnsemblAPI).get_info
                    if self.ensg_id == "" or self.ensg_id is None:
                        self.ensg_id = enst_dict.get("ensg_id")
                    if self.description == "" or self.description is None:
                        self.description = enst_dict.get("description")
                    if self.enst_id == "" or self.enst_id is None:
                        self.enst_id = enst_dict.get("enst_id")
                    if self.refseq_nt_id == "" or self.refseq_nt_id is None:
                        self.refseq_nt_id == enst_dict.get("refseq_nt_id")
                    if self.uniprot_id == "" or self.uniprot_id is None:
                        uniprot_id = enst_dict.get("uniprot_id")
                        if uniprot_id is not None and uniprot_id != "":
                            self.uniprot_id = enst_dict.get("uniprot_id")
            else:
                self.genename = human_ortholog_gene_id

        if self.genename != None:
            logger.info(f"Fetched orthologs for: {self.genename}")
        elif self.uniprot_id != None:
            logger.info(f"Fetched orthologs for: {self.uniprot_id}")
            
        self.had_orthologs_computed = True

    async def fetch_ortholog_async_semaphore(
        self,
        session: aiohttp.ClientSession,
        semaphore: asyncio.Semaphore,
        goaf: GOAnnotationsFile,
        human_ortholog_finder: Optional[HumanOrthologFinder] = None,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
    ) -> None:
        async with semaphore:
            await self.fetch_ortholog_async(
                session=session,
                goaf=goaf,
                human_ortholog_finder=human_ortholog_finder,
                uniprot_api=uniprot_api,
                ensembl_api=ensembl_api,
            )

    def fetch_info(
        self,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        required_keys=["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"],
    ) -> None:
        """
        Fetches additional information about this product. Additional information can be fetched if the Product has one of the four identifiers:
          - uniprot_id -> fetch info using uniprot_api.get_uniprot_info(self.uniprot_id) or ensembl_api.get_info(self.uniprot_id)
          - ensg_id -> fetch info using ensembl_api.get_info(self.ensg_id)
          - genename -> fetch info using ensembl_api.get_info(self.genename)

        The code checks the values for each Product's attribute from 'required_keys'. If any attributes are None, then
        the algorithm will attempt to find that information using queries in the following order:
          - uniprot_api.get_uniprot_info(self.uniprot_id) if uniprot_id != None
          - ensembl_api.get_info(self.ensg_id) if ensg_id != None
          - ensembl_api.get_info(self.genename) if genename != None
          - ensembl_api.get_info(self.uniprot_id) if uniprot_id != None

        After each query above, the returned dictionaries are processed and the attributes are set using
        setattr(self, key, value).

        Ideally, this function updates the following attributes: "genename", "description", "ensg_id", "enst_id", "refseq_nt_id"
        """
        self.had_fetch_info_computed = True
        if not (self.uniprot_id or self.genename or self.ensg_id):
            logger.debug(
                f"Product with id synonyms {self.id_synonyms} did not have an"
                " uniprot_id, gene name or ensg id. Aborting fetch info operation."
            )
            return
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        # required_keys = ["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"]
        # [TODO] Is uniprot really necessary? If it is faster, perhaps get uniprotID from genename and then first try to get info from uniprot
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.uniprot_id:
            info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
            if info_dict is not None:
                for key, value in info_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.ensg_id:
            enst_dict = ensembl_api.get_info(self.ensg_id)
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.genename:
            enst_dict = ensembl_api.get_info(self.genename)
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.uniprot_id:
            enst_dict = ensembl_api.get_info(self.uniprot_id)
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)

        # TODO: logger output which values are still missing

    async def fetch_info_async(
        self,
        client_session: aiohttp.ClientSession,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        required_keys=["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"],
    ) -> None:
        """
        required_keys correspond to the Product's attributes (class variables) that are checked. If any are None, then API requests
        are made so as to populate these variables with correct data.
        """
        self.had_fetch_info_computed = True
        if not (self.uniprot_id or self.genename or self.ensg_id):
            logger.debug(
                f"Product with id synonyms {self.id_synonyms} did not have an"
                " uniprot_id, gene name or ensg id. Aborting fetch info operation."
            )
            return
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.uniprot_id:
            info_dict = await uniprot_api.get_uniprot_info_async(
                self.uniprot_id, session=client_session
            )
            if info_dict is not None:
                for key, value in info_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.ensg_id:
            enst_dict = await ensembl_api.get_info_async(
                self.ensg_id, session=client_session
            )
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.genename:
            enst_dict = await ensembl_api.get_info_async(
                self.genename, session=client_session
            )
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)
        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.uniprot_id:
            enst_dict = await ensembl_api.get_info_async(
                self.uniprot_id, session=client_session
            )
            if enst_dict is not None:
                for key, value in enst_dict.items():
                    if value is not None and value != "" and value != []:
                        current_attr_value = getattr(self, key)
                        if (
                            current_attr_value is None
                            or current_attr_value == ""
                            or current_attr_value == []
                        ):
                            setattr(self, key, value)

    async def fetch_info_async_semaphore(
        self,
        session: aiohttp.ClientSession,
        semaphore: asyncio.Semaphore,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        required_keys=["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"],
    ):
        async with semaphore:
            await self.fetch_info_async(
                session, uniprot_api, ensembl_api, required_keys
            )

    def fetch_mRNA_sequence(self, ensembl_api: EnsemblApi) -> None:
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        sequence = ensembl_api.get_sequence(
            self.enst_id
        )  # enst_id because we want the mRNA transcript
        if sequence is not None:
            self.mRNA = sequence
        else:
            self.mRNA = -1

    @classmethod
    def from_dict(cls, d: dict) -> "Product":
        """
        Class method to create a new Product instance from a dictionary.

        Args:
            d (dict): The dictionary containing the data to create the Product instance.

        Returns:
            Product: A new Product instance created from the input dictionary.
        """
        return cls(
            d.get("id_synonyms"),
            d.get("genename"),
            d.get("uniprot_id"),
            d.get("description"),
            d.get("ensg_id"),
            d.get("enst_id"),
            d.get("refseq_nt_id"),
            d.get("mRNA"),
            d.get("scores") if "scores" in d else None,
            d.get("had_orthologs_computed") if "had_orthologs_computed" in d else False,
            (
                d.get("had_fetch_info_computed")
                if "had_fetch_info_computed" in d
                else False
            ),
        )

from typing import Optional, List, Union
import asyncio
import aiohttp
import copy

from .ModelSettings import ModelSettings
from .ModelStats import ModelStats
from ..parse.OrthologParsers import HumanOrthologFinder
from ..parse.GOAFParser import GOAnnotationsFile
from ..web_apis.EnsemblApi import EnsemblApi
from ..web_apis.UniprotApi import UniProtApi
from ..util.WebsiteParser import WebsiteParser
from ..util.ApiUtil import EnsemblUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class Product:
    def __init__(
        self,
        id_synonyms: List[str],
        taxon: str = None,
        target_taxon: str = None,
        genename: str = None,
        uniprot_id: str = None,
        description: str = None,
        ensg_id: Union[str,list] = None,
        enst_id: str = None,
        refseq_nt_id: str = None,
        mRNA: str = None,
        scores: dict = None,
        had_orthologs_computed: bool = False,
        had_fetch_info_computed: bool = False,
        gorth_ortholog_exists: bool = None,
        gorth_ortholog_status: str = None
    ):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id_synonyms (str): The list of ID of the product and synonyms. -> after ortholog translation it turns out that some products are the same. Example: RGD:3774, Xenbase:XB-GENE-5818802, UniProtKB:Q9NTG7
            taxon (str)
            target_taxon (str): The target taxon number for homology search
            genename (str)
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            ensg_id (str): Ensembl gene ID or gene IDs
            enst_id (str): Ensembl transcript ID.
            refseq_nt_id (str): Refseq (reference sequence) transcript ID.
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
            had_orthologs_computed (bool): If this Product instance has had the fetch_ortholog function called already.
            had_fetch_info_computed (bool): If this Product instance has had the fetch_info function called already.
            gorth_ortholog_exists (bool): True, if gOrth (from gProfiler) finds either 1 or multiple ortholog ids, otherwise false (if no orthologs are found via gOrth). None if gOrth query wasn't performed.
            gorth_ortholog_status (str): "definitive" for definitive gOrth orthologs (gOrth returns only one ortholog), "indefinitive" (gOrth returns multiple orthologs) or "none" if gOrth returns no orthologs
        """
        self.id_synonyms = id_synonyms
        self.taxon = taxon
        self.target_taxon = target_taxon
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
        self.gorth_ortholog_exists = gorth_ortholog_exists
        self.gorth_ortholog_status = gorth_ortholog_status

        # see if UniProtKB id is already in id_synonyms:
        for id_syn in self.id_synonyms:
            if "UniProt" in id_syn:
                self.uniprot_id = id_syn

    def update(self, other_product):
        """
        Updates this instance with new values.
          - id synonyms are appended
          - other values are updated if not None
        """
        assert(isinstance(other_product, Product))
        for attr_name in dir(self):
            if not callable(getattr(self, attr_name)) and not attr_name.startswith("__"):
                if attr_name == "id_synonyms": # extend and skip loop to prevent resetting previous values
                    self.id_synonyms.extend(other_product.id_synonyms)
                    continue
                if getattr(self, attr_name) is None:
                    setattr(self, attr_name, getattr(other_product, attr_name))

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
        # TODO: this function isn't updated to the latest pipeline; update to be the same as fetch_ortholog_async

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
                human_ortholog_gene_id = human_ortholog_finder.find_human_ortholog(self.id_synonyms[0])
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
                        info_dict = uniprot_api.get_uniprot_info(self.id_synonyms[0], taxon_id_num=self.taxon)  # bugfix
                    else:
                        info_dict = {"genename": goaf.get_uniprotkb_genename(self.id_synonyms[0])}
                else:  # self.uniprot_id exists
                    if _d_compare_goaf is True:
                        info_dict = uniprot_api.get_uniprot_info(self.uniprot_id, taxon_id_num=self.taxon)
                    else:
                        info_dict = {"genename": goaf.get_uniprotkb_genename(self.uniprot_id)}
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
                if offline_queried_ortholog is None:  # if algorithm enters this section due to _d_compare_goaf == True, then this accelerates code, as it prevents double calculations
                    human_ortholog_gene_id = human_ortholog_finder.find_human_ortholog(self.id_synonyms[0])  # file-based search; alternative spot for GOAF analysis
                else:
                    human_ortholog_gene_id = offline_queried_ortholog
                if human_ortholog_gene_id is None:  # if file-based search finds no ortholog
                    logger.warning(f"human ortholog finder did not find ortholog for {self.id_synonyms[0]}")
                    if self.ensg_id is not None:
                        if "ENSG" in self.ensg_id:
                            human_ortholog_gene_ensg_id = self.ensg_id
                        else:
                            human_ortholog_gene_ensg_id = ensembl_api.get_human_ortholog(self.ensg_id)
                    else:
                        human_ortholog_gene_ensg_id = ensembl_api.get_human_ortholog(self.id_synonyms[0])

                    if human_ortholog_gene_ensg_id is not None:
                        enst_dict = ensembl_api.get_info(human_ortholog_gene_ensg_id) # TODO: add taxon here!!
                        human_ortholog_gene_id = enst_dict.get("genename")
                        if human_ortholog_gene_id is not None:
                            if (
                                DROP_MIRNA_FROM_ENSEMBL_QUERY is True
                                and "MIR" in human_ortholog_gene_id
                            ):
                                human_ortholog_gene_id = None  # Ensembl query returned a miRNA, return None
        
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
                        if self.genename != human_ortholog_gene_id:  # with the current workflow, these will always be the same
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
        target_organism_id_number:Optional[int] = None,
        human_ortholog_finder: Optional[HumanOrthologFinder] = None,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
        model_settings: Optional[ModelSettings] = None
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
        initial_taxon = self.taxon
        initial_genename = self.genename
        initial_uniprot_id = self.uniprot_id
        initial_ensg_id = self.ensg_id
        starting_state = copy.copy(self)
        ortholog_query_pathway = ""
        comments = []

        if not human_ortholog_finder:
            human_ortholog_finder = HumanOrthologFinder()
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()
        
        if target_organism_id_number is None:
            if model_settings is None:
                if self.target_taxon is None:
                    raise Exception(f"target_organism_id_number or a ModelSettings instance was not passed to the ortholog query function. The target species cannot be assumed!")
        
        # determine target_organism_id_num
        target_organism_id_num = None
        if target_organism_id_number is not None:
            target_organism_id_num = target_organism_id_number
        if model_settings is not None:
            # check for match if model_settings and target_organism_id_number are both supplied as params
            if model_settings.target_organism is not None:
                if target_organism_id_num is not None:
                    if model_settings.target_organism.ncbi_id != target_organism_id_num:
                        logger.warning(f"ModelSettings target organism id number {model_settings.target_organism.ncbi_id} doesn't match with target_organism_id_num {target_organism_id_num} supplied to the ortholog query function. ModelSettings id number will be used for ortholog query.")
                target_organism_id_num = model_settings.target_organism.ncbi_id
        if self.target_taxon is not None and target_organism_id_num is None:
            if ":" in self.target_taxon:
                target_organism_id_num = self.target_taxon.split(":")[1]
            else:
                target_organism_id_num = self.target_taxon
        
        if any("UniProtKB" in syn for syn in self.id_synonyms) and f"{target_organism_id_num}" in self.taxon:
            # no ortholog query is performed on Homo Sapiens, just info fetch
            ortholog_query_pathway = "None"
            comments.append("UniProtKB matching target organism taxon in id synonyms - no ortholog query needed.")

            uniprot_id = next((syn for syn in self.id_synonyms if "UniProtKB" in syn), None) # select the UniProtKB id
            # fallback: check also self.uniprot_id
            if uniprot_id is None and self.uniprot_id is not None:
                uniprot_id = self.uniprot_id

            # perform online or offline query
            if(
                model_settings is not None
                and model_settings.uniprotkb_genename_online_query == True
            ):
                info_dict = await uniprot_api.get_uniprot_info_async(uniprot_id=uniprot_id, session=session, organism_taxon_id_num=self.taxon)
            else:
                info_dict = {"genename": goaf.get_uniprotkb_genename(uniprot_id)} # TODO: implement goafmaster -> human!!!
            self.genename = info_dict.get("genename")
            # TODO: query other data from info_dict, if it exists !
        else:
            # checked here are: UniProtKB not from the homo sapiens taxon, ZFIN, MGI, RGD, Xenbase, etc.
            # autoselect the first id synonym or uniprotkb, uniprotkb has precedence
            initial_gene_id = None
            ortholog_gene_ensg_id = None

            # UniProtKB ids don't work with ensembl_api.get_ortholog_async. If UniProtKB id is the sole id in id_synonyms, attempt
            # to find corresponding ENSG id
            if any("UniProtKB" in syn for syn in self.id_synonyms):
                if len(self.id_synonyms) == 1:
                    # UniProtKB is the only id synonym -> attempt ENSG
                    if self.ensg_id is not None:
                        initial_gene_id = self.ensg_id
                    else:
                        # attempt to convert uniprot id to ensembl
                        if self.gorth_ortholog_status != "none":
                            uniprot_lookup = await uniprot_api.get_uniprot_info_async(self.id_synonyms[0], session=session, organism_taxon_id_num=self.taxon)
                            if uniprot_lookup is None or uniprot_lookup == {}: # if conversion unsuccessful (and uniprotkb is the only id), dont query ortholog, because it will fail
                                return None
                            else: # if conversion was successful -> scoop data from dict
                                self.genename = uniprot_lookup.get("genename")
                                self.description = uniprot_lookup.get("description")
                                self.ensg_id = uniprot_lookup.get("ensg_id")
                                self.enst_id = uniprot_lookup.get("enst_id")
                                self.refseq_nt_id = uniprot_lookup.get("refseq_nt_id")
                                if self.ensg_id is None:
                                    # if ensg_id wasn't found (and uniprotkb is the only id), dont' query ortholog, because it will fail
                                    return None
                                else:
                                    # if ensg_id was found, use it to query ortholog
                                    initial_gene_id = self.ensg_id
                else:
                    # pick id synonym other than UniProtKB
                    for id_syn in self.id_synonyms:
                        if "UniProtKB" not in id_syn:
                            initial_gene_id = id_syn
                            break
            else:
                # autopick the first id synonym
                initial_gene_id = self.id_synonyms[0]
            
            # check gOrth ortholog query status (gOrth batch ortholog query can be performed prior to a regular ortholog query)
            # if gOrth didn't discover ortholog, and model_settings.gorth_ortholog_refetch is False, then don't query ortholog
            should_query_ortholog = True
            if (
                model_settings is not None
                and model_settings.gorth_ortholog_refetch == False
                and self.gorth_ortholog_exists == False
            ): # if gOrth ortholog query did not find orthologs, but gorth_ortholog_refetch is False in model settings -> don't fetch
                should_query_ortholog = False
                ortholog_query_pathway = "None"
                comments.append("gOrth gene ortholog query didn't find ortholog, model settings gorth_ortholog_refetch is False.")
                # return None
            
            # check ensg id (this can be computed by gOrth batch ortholog query)
            if (
                self.ensg_id is not None 
                and self.gorth_ortholog_exists == True
            ):
                # this is either a definitive ortholog queried by gOrth or an indefinitive ortholog queried by gOrth using model_settings.gorth_ortholog_fetch_for_indefinitive_orthologs = False -> no need to query ortholog
                # if model_settings.gorth_ortholog_fetch_for_indefinitive_orthologs is True, then
                # self.ensg_id will NOT be set (in the gOrth query function), thus ortholog will be queried in the regular pipeline.
                ortholog_query_pathway = "Pre-existing gOrth query"
                comments.append("gOrth already found ortholog for this gene.")
                should_query_ortholog = False
                ortholog_gene_ensg_id = self.ensg_id
            
            # get the stable id prefix of the target organism. If target organism is zebrafish, then stable id prefix is "ENSDAR".
            # If the Ensembl id of this product is already the same as stable id prefix (in our case, ENSDAR), then don't query ortholog.

            target_organism_stable_id_prefix = WebsiteParser.get_ensembl_stable_id_prefixes_table()[f"{target_organism_id_num}"]['stable_id_prefix'] # ex. "ENSDAR" for danio rerio, "ENS" for human
            if self.ensg_id is not None:
                current_organism_stable_id_prefix = EnsemblUtil.split_ensembl_id(self.ensg_id)['stable_id_prefix']
                stable_id_prefix_match = (target_organism_stable_id_prefix == current_organism_stable_id_prefix)
                if stable_id_prefix_match:
                    # (Product).fetch_ortholog_async was called with same target organism as is the organism assigned to this product -> this ortholog has already been found!
                    ortholog_query_pathway = "None"
                    comments.append(f"Existing ensg id stable prefix {current_organism_stable_id_prefix} is the same as target organism stable id prefix.")
                    should_query_ortholog = False
                    ortholog_gene_ensg_id = self.ensg_id
            
            if should_query_ortholog:
                # attempt offline, file-based search
                human_ortholog_gene_symbol = None
                if self.id_synonyms[0].split(":")[0] in HumanOrthologFinder.get_supported_organism_dbs():
                    ortholog_query_pathway = "HumanOrthologFinder"
                    human_ortholog_gene_symbol = (await human_ortholog_finder.find_human_ortholog_async(self.id_synonyms[0]))

                # attempt online Ensembl query
                if human_ortholog_gene_symbol is None:
                    logger.info(f"HumanOrthologFinder did not find gene symbol for {self.id_synonyms[0]}. Attempting Ensembl ortholog query.")
                    # ortholog_gene_ensg_id = await ensembl_api.get_human_ortholog_async(id=initial_gene_id, session=session, taxon=f"NCBITaxon:{target_organism_id_num}") # replaced to offer full taxon modularity
                    ortholog_gene_ensg_id = await(ensembl_api.get_ortholog_async(id=initial_gene_id, session=session, source_taxon=self.taxon, target_taxon=f"NCBITaxon:{target_organism_id_num}"))
                    ortholog_query_pathway = f"EnsemblApi (get_human_ortholog: id = {initial_gene_id}, taxon = NCBITaxon:{target_organism_id_num})"
            else:
                # ensg id already exists
                if self.ensg_id is not None:
                    ortholog_gene_ensg_id = self.ensg_id
            
            # ortholog_gene_ensg_id should now be filled out
            if ortholog_gene_ensg_id is not None:
                current_stable_id_prefix = EnsemblUtil.split_ensembl_id(ortholog_gene_ensg_id)['stable_id_prefix']
                stable_id_prefix_match = (current_stable_id_prefix == target_organism_stable_id_prefix)

                if stable_id_prefix_match == True:
                    ortholog_query_pathway += f" + EnsemblApi (get_info: {ortholog_gene_ensg_id})"
                    enst_dict = await ensembl_api.get_info_async(ortholog_gene_ensg_id, taxon=self.taxon, session=session)
                    self.genename = enst_dict.get("genename")

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
                    # self.genename = ortholog_gene_ensg_id
                    pass

        ModelStats.product_ortholog_query_results[self.id_synonyms[0]] = {
            'initial_state': starting_state.to_json(),
            'end_state': self.to_json(),
            'delta_genename': self.genename if initial_genename != self.genename else False,
            'delta_uniprot_id': self.uniprot_id if initial_uniprot_id != self.uniprot_id else False,
            'delta_ensg_id': self.ensg_id if initial_ensg_id != self.ensg_id else False,
            'ortholog_query_pathway': ortholog_query_pathway,
            'comments': comments
        }

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
        target_organism_taxon_number:int=9606,
        model_settings: Optional[ModelSettings] = None,
        human_ortholog_finder: Optional[HumanOrthologFinder] = None,
        uniprot_api: Optional[UniProtApi] = None,
        ensembl_api: Optional[EnsemblApi] = None,
    ) -> None:
        async with semaphore:
            await self.fetch_ortholog_async(
                session=session,
                goaf=goaf,
                target_organism_id_number=target_organism_taxon_number,
                model_settings=model_settings,
                human_ortholog_finder=human_ortholog_finder,
                uniprot_api=uniprot_api,
                ensembl_api=ensembl_api
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
            info_dict = uniprot_api.get_uniprot_info(self.uniprot_id, taxon_id_num=self.taxon)
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
            enst_dict = ensembl_api.get_info(self.ensg_id, taxon=self.target_taxon) # ENSG id should have been converted to target_taxon by now
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
            enst_dict = ensembl_api.get_info(self.genename, taxon=self.taxon) # for genename, attempt source taxon!
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
            enst_dict = ensembl_api.get_info(self.uniprot_id, taxon=self.taxon)
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
            logger.debug(f"Product with id synonyms {self.id_synonyms} did not have an uniprot_id, gene name or ensg id. Aborting fetch info operation.")
            return
        if not uniprot_api:
            uniprot_api = UniProtApi()
        if not ensembl_api:
            ensembl_api = EnsemblApi()

        if (
            any(getattr(self, key) is None for key in required_keys)
            or any(getattr(self, key) == "" for key in required_keys)
        ) and self.uniprot_id:
            info_dict = await uniprot_api.get_uniprot_info_async(self.uniprot_id, session=client_session, organism_taxon_id_num=self.taxon)
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
            enst_dict = await ensembl_api.get_info_async(self.ensg_id, taxon=self.taxon, session=client_session)
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
            enst_dict = await ensembl_api.get_info_async(self.genename, taxon=self.taxon, session=client_session)
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
            enst_dict = await ensembl_api.get_info_async(self.uniprot_id, taxon=self.taxon, session=client_session)
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
            d.get("taxon") if "taxon" in d else None,
            d.get("target_taxon") if "target_taxon" in d else None,
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
    
    def to_json(self):
        """
        Constructs a JSON representation of this class. This is used during (ReverseLookup).save_model to save the ModelSettings
        """
        json_data = {}
        for attr_name, attr_value in vars(self).items():
            # custom handling for target_organism and ortholog_organisms, as they are code objects -> convert them to json
            if not callable(attr_value) and not attr_name.startswith("__"):
                # append to json_data result dict
                json_data[attr_name] = attr_value
        return json_data

from typing import List
import os
import gzip
import urllib

from ..core.ModelSettings import ModelSettings
from .OBOParser import OboParser
from ..util.FileUtil import FileUtil
from ..util.DictUtil import DictUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger

class GOAFMaster:
    def __init__(
            self,
            goa_filepaths: dict = {
                'homo_sapiens': {
                    'organism': "homo_sapiens",
                    'ncbi_taxon_id': "9606",
                    'local_filepath': "data_files/goa_human.gaf",
                    'download_url': "http://geneontology.org/gene-associations/goa_human.gaf.gz",
                },
                'danio_rerio': {
                    'organism': "danio_rerio",
                    'ncbi_taxon_id': "7955",
                    'local_filepath': "data_files/zfin.gaf",
                    'download_url': "http://current.geneontology.org/annotations/zfin.gaf.gz"
                },
                'rattus_norvegicus': {
                    'organism': "rattus_norvegicus",
                    'ncbi_taxon_id': "10116",
                    'local_filepath': "data_files/rgd.gaf",
                    'download_url': "http://current.geneontology.org/annotations/rgd.gaf.gz"
                },
                'mus_musculus': {
                    'organism': "mus_musculus",
                    'ncbi_taxon_id': "10090",
                    'local_filepath': "data_files/mgi.gaf",
                    'download_url': "http://current.geneontology.org/annotations/mgi.gaf.gz"
                },
                'xenopus': {
                    'organism': "xenopus",
                    'ncbi_taxon_id': "8353",
                    'local_filepath': "data_files/xenbase.gaf",
                    'download_url': "http://current.geneontology.org/annotations/xenbase.gaf.gz"
                }
            },
            go_categories: list = ["biological_process", "molecular_activity", "cellular_component"],
            valid_evidence_codes:list = [],
            evidence_codes_to_ecoids:dict = {}
    ):
        """
        This is a master-class, which allows managing of several GOAnnotationsFiles (for different organisms). Parameters:
          - goa_filepaths: a dictionary of the supplied GOA files. In it's minimum, this dictionary has to be structured like the following:

            {
                "GOA_FILE_NAME": {
                    'organism': "ORGANISM_NAME",
                    'ncbi_taxon_id': "NCBI_TAXON_ID",
                    'local_filepath': "LOCAL_FILEPATH",
                    'download_url': "DOWNLOAD_URL"
                }
            }

            Example:

            {
                "homo_sapiens": {
                    'organism': "homo_sapiens",
                    'ncbi_taxon_id': 9606
                    'local_filepath': "data_files/goa_human.gaf",
                    'download_url': "http://geneontology.org/gene-associations/goa_human.gaf.gz"
                }
            }

            If GOAFMaster were supplied with the above 'goa_filepaths' parameter, it would append the created GOA instance to self.goa_files,
            which could be accessed by self.goa_files[ORGANISM], in our case self.goa_files["homo_sapiens"].

          - go_categories: a list of GO categories to take into account. They must be supplied in the list as strings. Available options are: "biological_process", "molecular_activity" and "cellular_component"
          - (list) valid_evidence_codes: a list of valid evidence codes. Annotations which do not have an evidence code aligned with this list will be discarded.
                                         WARNING: supplied evidence codes must be in the user-friendly format (e.g. IEA, TAS, ...) and not the ECO format! If they are supplied in the ECO format, a conversion attempt will be made using 'evidence_codes_to_ecoids' parameter, if set
          - (dict) evidence_codes_to_ecoids: a mapping from user-friendly evidence codes to their respective ECO ids. It is advisable to set this as a fallback if supplied evidence codes are in the ECO format.

        """
        self.goa_filepaths = goa_filepaths
        self.go_categories = go_categories
        self.goa_files = {} # dict mapping organism labels to respective goa files

        for organism_key, file_raw_info in goa_filepaths.items():
            goa_file = GOAnnotationsFile(
                filepath=file_raw_info['local_filepath'],
                download_url=file_raw_info['download_url'],
                go_categories=go_categories,
                organism_label=file_raw_info['organism'],
                organism_ncbi_taxon_id=file_raw_info['ncbi_taxon_id'],
                valid_evidence_codes=valid_evidence_codes,
                evidence_codes_to_ecoids=evidence_codes_to_ecoids
            )
            self.goa_files[organism_key] = goa_file

    def get_all_products_for_goterm(self, goterm_id:str, organism_labels=['homo_sapiens', 'danio_rerio', 'rattus_norvegicus', 'mus_musculus', 'xenopus'], indirect_annotations:bool = False, obo_parser:OboParser=None):
        """
        Finds all genes associated to goterm_id across differeng GAF files (specified by organims_labels, which should match the organism labels used to create the GOAFMaster instance).
        
        Returns:
          - [0]: a list of all gene ids associated with goterm_id from all of the specified organism_labels (all of the specified gaf files)
          - [1]: a dict mapping organism label keys to found gene ids

        Example usage:
        (GOAFMaster).get_all_products_for_goterm(GOID, ["homo_sapiens", "danio_rerio"])
        -> returns: 
             [0]: [ABC1, VEGFA, ..., ijp2, tsp5, mmp1, ...]
             [1]: {
                    'homo_sapiens': [ABC1, VEGFA],
                    'danio_rerio': [ijp2, tsp5, mmp1]
                  }
        """
        annotated_genes_simple_list = [] # list of annotated genes to goterm_id; source organism isn't known
        annotated_genes_dict = {} # dict between organism label (key) and associated genes
        for organism_label in organism_labels:
            if organism_label not in self.goa_files:
                logger.warning(f"Organism label '{organism_label}' was NOT found in self.goa_files. Was the correct GOA file created during construction of GOAFMaster?")
                continue
            goa_current = self.goa_files[organism_label]
            assert isinstance(goa_current, GOAnnotationsFile)
            current_genes = goa_current.get_all_products_for_goterm(goterm_id=goterm_id, indirect_annotations=indirect_annotations, obo_parser=obo_parser)
            annotated_genes_simple_list = [*annotated_genes_simple_list, *current_genes]
            annotated_genes_dict[organism_label] = current_genes

        return [annotated_genes_simple_list, annotated_genes_dict]
    
    def get_all_goterms_for_product(self, product_id: str, organism_labels=['homo_sapiens', 'danio_rerio', 'rattus_norvegicus', 'mus_musculus', 'xenopus'], indirect_annotations:bool=False, obo_parser:OboParser=None):
        """
        Finds all GO Terms associated to product_id across differeng GAF files (specified by organims_labels, which should match the organism labels used to create the GOAFMaster instance).
        
        Returns:
          - [0]: a list of all GO Term ids associated with product_id from all of the specified organism_labels (all of the specified gaf files)
          - [1]: a dict mapping organism label keys to found GO Term ids

        Example usage:
        (GOAFMaster).get_all_goterms_for_product(gene_id, ["homo_sapiens", "danio_rerio"])
        -> returns: 
             [0]: [ABC1, VEGFA, ..., ijp2, tsp5, mmp1, ...]
             [1]: {
                    'homo_sapiens': [ABC1, VEGFA],
                    'danio_rerio': [ijp2, tsp5, mmp1]
                  }
        """
        annotated_termids_simple_list = []
        annotated_termids_dict = {}
        
        for organism_label in organism_labels:
            if organism_label not in self.goa_files:
                logger.warning(f"Organism label '{organism_label}' was NOT found in self.goa_files. Was the correct GOA file created during construction of GOAFMaster?")
                continue
            goa_current = self.goa_files[organism_label]
            assert isinstance(goa_current, GOAnnotationsFile)
            current_termids = goa_current.get_all_terms_for_product(product=product_id, indirect_annotations=indirect_annotations, obo_parser=obo_parser)
            annotated_termids_simple_list = [*annotated_termids_simple_list, *current_termids]
            annotated_termids_dict[organism_label] = current_termids

        return [annotated_termids_simple_list, annotated_termids_dict]
    
class GOAnnotationsFile:
    def __init__(
        self,
        filepath: str = "data_files/goa_human.gaf",
        download_url = "http://geneontology.org/gene-associations/goa_human.gaf.gz",
        go_categories: list = [
            "biological_process",
            "molecular_activity",
            "cellular_component",
        ],
        organism_label = "",
        organism_ncbi_taxon_id = "",
        valid_evidence_codes:list = [],
        evidence_codes_to_ecoids:dict = {}
    ) -> None:
        """
        TODO: update comment

        This class provides access to a Gene Ontology Annotations File, which stores the relations between each GO Term and it's products (genes),
        along with an evidence code, confirming the truth of the interaction. A GO Annotation comprises of a) GO Term, b) gene / gene product c) evidence code.
        WARNING: GOAF stores only DIRECT annotations (See https://geneontology.org/docs/faq/ "Why does AmiGO display annotations to term X but these annotations aren’t in the GAF file?")

        Parameters:
          - (str) filepath: the filepath to the GO Annotations File downloaded file from http://current.geneontology.org/products/pages/downloads.html -> Homo Sapiens (EBI Gene Ontology Database) - protein = goa_human.gaf; link = http://geneontology.org/gene-associations/goa_human.gaf.gz
                            if left to default value, self._filepath will be set to 'app/goreverselookup/data_files/goa_human.gaf'. The file should reside in app/goreverselookup/data_files/ and the parameter filepath should be the file name of the downloaded file inside data_files/
          - (str) download_url: the download url that is used to download this file from the web
          - (list) go_categories: determines which GO categories are valid. Default is that all three GO categories are valid. Setting GO categories determines which products
                                  or terms are returned from goaf.get_all_products_for_goterm and goaf.get_all_terms_for_product functions. The algorithm excludes any associations whose category doesn't match go_categories already in
                                  the GOAF file read phase - lines not containing a desired category (from go_categories) won't be read.
          - (str) organism_label: a descriptive label of the organism, for which this GAF is constructed (eg. homo_sapiens)
          - (str) organism_ncbi_taxon_id: ncbi taxon id for the organism, for which this GAF is constructed(eg. 9606 for homo_sapiens)
          - (list) valid_evidence_codes: a list of valid evidence codes. Annotations which do not have an evidence code aligned with this list will be discarded.
                                         WARNING: supplied evidence codes must be in the user-friendly format (e.g. IEA, TAS, ...) and not the ECO format! If they are supplied in the ECO format, a conversion attempt will be made using 'evidence_codes_to_ecoids' parameter, if set
          - (dict) evidence_codes_to_ecoids: a mapping from user-friendly evidence codes to their respective ECO ids. It is advisable to set this as a fallback if supplied evidence codes are in the ECO format.
                                  
        See also:
          - http://geneontology.org/docs/download-go-annotations/
          - http://current.geneontology.org/products/pages/downloads.html
        """
        # check that all valid_evidence_codes are in the user-friendly (and not in the ECO) format
        if evidence_codes_to_ecoids != {}:
            _v = []
            ecoids_to_evidence_codes = DictUtil.reverse_dict(evidence_codes_to_ecoids)
            for vec in valid_evidence_codes:
                if "ECO:" in vec:
                    converted = ecoids_to_evidence_codes.get(vec, None)
                    if vec is not None:
                        _v.append(converted)
                else:
                    _v.append(vec)
            logger.debug(f"Converted evidence codes:")
            logger.debug(f"  - input: {valid_evidence_codes}")
            logger.debug(f"  - output: {_v}")
            valid_evidence_codes = _v

        self.go_categories = go_categories
        self.organism_label = organism_label
        self.organism_ncbi_taxon_id = organism_ncbi_taxon_id

        if filepath == "" or filepath is None:
            self._filepath = "data_files/goa_human.gaf"
        else:
            self._filepath = filepath

        logger.info(f"Attempting to create GOAF using: {self._filepath}")

        # self._check_file()
        if not os.path.exists(self._filepath) or FileUtil.is_file_empty(self._filepath):
            FileUtil.download_zip_file(filepath=self._filepath, download_url=download_url, zip_specifier="rt")

        logger.info("  - GOAF filepath exists.")
        with open(self._filepath, "r") as read_content:
            temp_content = read_content.readlines()
            self._readlines = []
            for line in temp_content:
                if not line.startswith("!") and not line.strip() == "":
                    line = line.strip()
                    line_category = self._get_go_category_from_line(line)
                    line_evidence_code = self._get_evidence_code_from_line(line)
                    if line_category in go_categories and line_evidence_code in valid_evidence_codes:
                        self._readlines.append(line)
        self.terms_dict = None
        self.products_dict = None
        logger.info(f"  - GOAF created with {len(self._readlines)} annotations.")

    # this is obsolete. TODO: remove after ensuring it doesn't cause problems
    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if os.path.exists(self._filepath):
            return True
        else:
            url = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
            # download the gzip file and save it to a temporary file
            temp_file, _ = urllib.request.urlretrieve(url)

            # read the contents of the gzip file and save it to the txt file
            with gzip.open(temp_file, "rt") as f_in, open(self._filepath, "w") as f_out:
                for line in f_in:
                    f_out.write(line)

            # delete the temporary file
            os.remove(temp_file)

        if os.path.exists(self._filepath):
            return True
        else:
            return False

    def _get_go_category_from_line(self, line: str):
        """
        Expects a line similar to:
            UniProtKB	A0A075B6H8	IGKV1D-42	involved_in	GO:0002250	GO_REF:0000043	IEA	UniProtKB-KW:KW-1064	P	Probable non-functional immunoglobulin kappa variable 1D-42	IGKV1D-42	protein	taxon:9606	20230306	UniProt
        and returns the GO aspect (GO category) of a line, which can be either:
            - P --> biological_process
            - F --> molecular_activity
            - C --> cellular_component

        A word about the search algorithm:
        Line elements in the GOAF are in the following order:
          - [0] DB
          - [1] DB Object Id
          - [2] Db Object Symbol
          - [3] Qualifier (optional)
          - [4] GO Id
          - [5] DB:Reference
          - [6] Evidence code
          - [7] With (or) From (optional)
          - [8] Aspect
          - [9] DB Object Name
          - [10] DB Object Synonym
          - [11] DB Object Type
          - [12] Taxon
          - [13] Date
          - [14] Assigned by
          - [15] Annotation extension
          - [16] Gene product form id
        Since two line elements (qualifier) and (with or from) are optional (before the Aspect element), the array element at index 8
        will be checked if it contains only one letter. If not, elements at index 7 and 6 will be checked if they contain only one letter (corresponding to Aspect),
        since one of these elements must have the Aspect. In other words, Aspect can be at indices:
          - 6: if both "qualifier" and "with or from" elements are missing
          - 7: if only one ("qualifier" or "with or from") is missing
          - 8: if "qualifier" and "with or from" elements are supplied.
        """
        line_split = []
        if isinstance(line, list):
            line_split = line
        else:
            # line is not split, split on tab
            line_split = line.split("\t")

        # aspect is the same as go_category
        aspect = ""
        start_index = 8  # expected Aspect index if line
        for i in range(3):  # will go: 0,1,2
            aspect_index = start_index - i  # will go: 8,7,6
            aspect = line_split[
                aspect_index
            ]  # query line elements 8, 7 and 6 (possible line locations for Aspect)
            if (
                len(aspect) == 1
            ):  # if length of Aspect string is 1, then Aspect is truly P, C or F
                break

        match aspect:
            case "P":
                return "biological_process"
            case "F":  # molecular function in https requests when querying GO Terms associated with gene ids is returned as molecular_activity
                return "molecular_activity"
            case "C":
                return "cellular_component"
        return None

    def _get_evidence_code_from_line(self, line):
        """
        Extracts evidence code from the line. It is at position 5 or 6 in the split line and is of length 2 or 3 characters.
        Evidence codes in the GOA files are in the user-friendly format (e.g. IEA) and not in the ECO identifier format.
        The element before evidence code is always a DB:reference (eg GO_REF:xxxx or Reactome:xxx, ...)

        Params:
          - line: either a line already split at (\t characters). If line is a str, it will be split at tabs.
        
        Returns: the evidence code of this line
        """
        if isinstance(line, str):
            line = line.split("\t")
        evidence_code = line[6]
        if len(evidence_code) == 2 or len(evidence_code) == 3:
            return evidence_code
        evidence_code = line[5]
        if len(evidence_code) == 2 or len(evidence_code) == 3:
            return evidence_code
        return None

    def get_all_products_for_goterm(self, goterm_id: str, model_settings:ModelSettings=None,indirect_annotations:bool=False, obo_parser:OboParser=None) -> List[str]:
        """
        This method returns all unique products associated with the GO term id.
        The return of this function is influenced by the go_categories supplied to the constructor of the GOAF!

        Args:
          - (str) goterm_id: a GO Term identifier, eg. GO:0003723
          - (ModelSettings) model_settings: must be passed if using indirect annotations to determine indirect annotations direction (ie children / parents) and the max indirect annotations depth
          - (bool) indirect_annotations: if True, will return a list of unique products for the specified goterm_id and for all of the indirect goterms of the 'goterm_id'
          - (OboParser) obo_parser: an OboParser instance, to prevent recalculations. If no OboParser is passed, then this function will attempt to call OboParser() to create an OboParser instance.

        Returns:
            List[str]: a List of all product (gene/gene products) gene names, eg. ['NUDT4B', ...]

        Example: for 'GO:0003723' it returns ['KMT2C', 'CLNS1A', 'ZCCHC9', 'TERT', 'BUD23', 'DDX19B', 'CCAR2', 'NAP1L4', 'SAMSN1', 'ERVK-9', 'COA6', 'RTF1', 'AHCYL1', 'SMARCA4', ... (total len = 1378)]
        """

        if self.terms_dict is None:
            self.populate_terms_dict()

        direct_products = self.terms_dict.get(goterm_id, [])
        if indirect_annotations == False:
            return direct_products
        else:
            if obo_parser is None:
                obo_parser = OboParser()

            indirect_products = []
            indirect_ids = obo_parser.get_indirect_annotations(
                term_id=goterm_id,
                indirect_annotations_direction=model_settings.indirect_annotations_direction,
                max_depth=model_settings.indirect_annotations_max_depth             
            )
            for indirect_id in indirect_ids:
                indirect_products = self.terms_dict.get(indirect_id, [])
                indirect_products += indirect_products
            
            return set(direct_products + indirect_products)
            

    def populate_poducts_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a product gene name and an associated GO Term.

        The result is a dictionary (self.products_dict), mapping keys (product gene names, eg. NUDT4B) to a List of all associated
        GO Terms (eg. ['GO:0003723', ...])
        """
        self.products_dict = {}
        for line in self._readlines:  # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split("\t")
            self.products_dict.setdefault(chunks[2], set()).add(chunks[4])  # create a key with the line's product gene name (if the key already exists, don't re-create the key - specified by the setdefault method) and add the associated GO Term to the value set. eg. {'NUDT4B': {'GO:0003723'}}, after first line is processed, {'NUDT4B': {'GO:0003723'}, 'NUDT4B': {'GO:0046872'}} after second line ...
        for key,values in self.products_dict.items():  # the set() above prevents the value elements (GO Terms) in dictionary to be repeated
            self.products_dict[key] = list(values)  # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}

    def populate_terms_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a GO Term and it's associated product gene name.

        The result is a dictionary (self.terms_dict), mapping keys (GO Terms, eg. GO:0003723) to a List of all
        associated product gene names (eg. ['NUDT4B', ...])
        """
        self.terms_dict = {}
        for line in self._readlines:  # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split("\t")
            self.terms_dict.setdefault(chunks[4], set()).add(chunks[2])  # create a key with the line's GO Term (if the key already exists, don't re-create the key - specified by the setdefault method) and add the product' gene name to the value set. eg. {'GO:0003723': {'NUDT4B'}}, after first line is processed, {'GO:0003723': {'NUDT4B'}, 'GO:0046872': {'NUDT4B'}} after second line ...
        for key,values in self.terms_dict.items():  # the previous set() prevents the value elements (product gene names) in dictionary to be repeated
            self.terms_dict[key] = list(values)  # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}

    def get_all_terms_for_product(self, product: str, model_settings:ModelSettings=None, indirect_annotations:bool=False, obo_parser:OboParser=None) -> List[str]:
        """
        Gets all GO Terms associated to a product gene name.
        The return of this function is influenced by the go_categories supplied to the constructor of the GOAF!

        Args:
          - (str) product: must be a gene name corresponding to a specific gene/gene product, eg. NUDT4B
          - (ModelSettings) model_settings: must be passed if using indirect annotations to determine indirect annotations direction (ie children / parents) and the max indirect annotations depth
          - (bool) indirect_annotations: if True, will also return all indirectly annotated terms for product (ie. all sub-terms of directly annotated terms in the GOAF)
          - (OboParser) obo_parser: an OboParser instance, to prevent recalculations. If no OboParser is passed, then this function will attempt to call OboParser() to create an OboParser instance.

        Returns:
          - List[str]: a List of all GO Term ids associated with the input product's gene name

        Example: for 'NUDT4B', it returns ['GO:1901911', 'GO:0071543', 'GO:0005737', 'GO:0000298', 'GO:0005634', 'GO:0034431', 'GO:0034432', 'GO:0046872', 'GO:0008486', 'GO:1901909', 'GO:0003723', 'GO:1901907', 'GO:0005829']
        """
        if self.products_dict is None:
            self.populate_poducts_dict()

        direct_annotations = self.products_dict.get(product, [])
        if indirect_annotations == False:
            return direct_annotations
        else:
            if obo_parser is None:
                obo_parser = OboParser()

            indirect_annotations = []
            for goterm_id in direct_annotations:
                indirect_annotations += obo_parser.get_indirect_annotations(
                    term_id=goterm_id,
                    indirect_annotations_direction=model_settings.indirect_annotations_direction,
                    max_depth=model_settings.indirect_annotations_max_depth
                )
            
            return (direct_annotations + indirect_annotations)


    def get_all_terms(self) -> List[str]:
        """
        Returns a List of all unique GO Terms read from the GO Annotations file.
        In the current (27_04_2023) GO Annotation File, there are 18880 unique GO Terms.

        The return of this function is influenced by the go_categories supplied to the constructor of the GOAF!
        """
        if not self.terms_dict:
            self.populate_terms_dict()

        terms_list = [k for k, v in self.terms_dict.items()]
        return terms_list

    """ # THIS IS INVALID, orthologs cannot be queried from the GOAF !!!
    def get_all_product_orthologs(self, product_id:str):
        #""
        #Gets all orthologs in line for a specific product (gene) id. This function uses GOAF for the ortholog query.
        #TODO: check if this is actually even valid.
        #""
        # if a user sends a uniprotkb product_id here, default to get_uniprotkb_genename
        if "UniProtKB" in product_id:
            genename =  self.get_uniprotkb_genename(product_id)
            if genename != None:
                return genename
        
        possible_orthologs = {} # a dict between possible orthologs and the readlines where they are found
        for line in self._readlines:
            if product_id in line:
                line_elements = line.split("\t")
                if line_elements[1] != product_id:
                    # query the 8th line element With (or) From
                    # GOAF line elements: (1) DB (2) DB Object Id (3) DB Object Symbol (4) Qualifier (optional) (5) GO ID (6) DB:Reference (7) Evidence Code (8) With (or) from (optional) (9) Aspect ...
                    # We are interested in the 8th line element, but it is optional. Furthermore, Qualifier (4th line element) is also optional.
                    # However, we are certain that the "With or from" line element will appear after the Evidence Code (which is always a 3-character code - http://geneontology.org/docs/guide-go-evidence-codes/) and before the
                    # Aspect, which can be either "P" (biological process), "F" (biological function) or "C" (cellular component). If the difference between the evidence code index and the aspect index (index = position in the array) is
                    # greater than 1, then we are sure that the element between them is a "With or from" element.
                    evidence_code = ""
                    aspect = ""
                    evidence_code_index = -1
                    aspect_index = -1
                    i=0
                    for line_element in line_elements:
                        if len(line_element) == 3: # this is the Evidence Code
                            evidence_code = line_element
                            evidence_code_index = i
                        if len(line_element) == 1 and i > evidence_code_index: # this is the Aspect
                            aspect = line_element
                            aspect_index = i
                        i+=1
                    if aspect_index - evidence_code_index > 1:
                        # orthologs exist
                        orthologs = 
        """

    def get_uniprotkb_genename(self, product_id: str):
        """
        Gets the gene name (DB Object Symbol) for the supplied UniProtKB product_id.

        Parameters:
          - product_id: must be in the format UniProtKB:XXXXX

        Algorithm:
            If the product_id is a UniProtKB, then the GOAF is browsed to obtain the gene name, which
            is the third line element (DB Object Symbol) in the GOAF. If the parse doesnt find the uniprot id
            (XXXXX) in the GOAF as the second line element, then the supplied UniProtKB may be an animal protein.
            In this case, the program also records any lines that don't have (XXXXX) as the second line element, but still
            contain the (XXXXX) in the line. The program reverts to these lines and attempts to find a human ortholog. If
            all lines result in the same human ortholog, then the search was successful. TODO: implement logic if multiple different
            orthologs are found.
        """
        if "UniProtKB" in product_id:
            product_id = product_id.split(":")[
                1
            ]  # gets the raw id; UniProtKB:XXXXX -> XXXXX
        else:
            logger.warning(
                f"get_uniprotkb_genename unsucessful for {product_id}. Product id must"
                " be supplied in the UniProtKB:XXXXX format!"
            )
            return None

        gene_name = ""
        ortholog_lines = (
            []
        )  # lines which contain product_id, but not as the second element
        for line in self._readlines:
            if product_id in line:
                line_elements = line.split("\t")
                if line_elements[1] == product_id:
                    gene_name = line_elements[2]
                    return gene_name
                elif line_elements[1] != product_id:
                    ortholog_lines.append(line)

        if gene_name == "" and len(ortholog_lines) > 0:
            # goaf file was read, but product_id was never the second line element,
            # but was found in some lines to be an ortholog to some genes? or is maybe involved in some genes?
            # TODO: implement logic

            # for ortho_line in ortholog_lines:
            #   ...

            gene_name = ""  # delete this

        if gene_name == "":
            # if gene_name is still not found, return None

            return None

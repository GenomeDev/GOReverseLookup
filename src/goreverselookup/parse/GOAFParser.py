from typing import List
import os
import gzip
import urllib

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class GOAnnotationsFile:
    def __init__(
        self,
        filepath: str = "app/goreverselookup/data_files/goa_human.gaf",
        go_categories: list = [
            "biological_process",
            "molecular_activity",
            "cellular_component",
        ],
    ) -> None:
        """
        This class provides access to a Gene Ontology Annotations File, which stores the relations between each GO Term and it's products (genes),
        along with an evidence code, confirming the truth of the interaction. A GO Annotation comprises of a) GO Term, b) gene / gene product c) evidence code.

        Parameters:
          - (str) filepath: the filepath to the GO Annotations File downloaded file from http://current.geneontology.org/products/pages/downloads.html -> Homo Sapiens (EBI Gene Ontology Database) - protein = goa_human.gaf; link = http://geneontology.org/gene-associations/goa_human.gaf.gz
                            if left to default value, self._filepath will be set to 'app/goreverselookup/data_files/goa_human.gaf'. The file should reside in app/goreverselookup/data_files/ and the parameter filepath should be the file name of the downloaded file inside data_files/
          - (list) go_categories: determines which GO categories are valid. Default is that all three GO categories are valid. Setting GO categories determines which products
                                  or terms are returned from goaf.get_all_products_for_goterm and goaf.get_all_terms_for_product functions. The algorithm excludes any associations whose category doesn't match go_categories already in
                                  the GOAF file read phase - lines not containing a desired category (from go_categories) won't be read.

        See also:
          - http://geneontology.org/docs/download-go-annotations/
          - http://current.geneontology.org/products/pages/downloads.html
        """
        self.go_categories = go_categories
        if filepath == "" or filepath is None:
            self._filepath = "app/goreverselookup/data_files/goa_human.gaf"
        else:
            self._filepath = filepath

        logger.info(f"Attempting to create GOAF using: {self._filepath}")

        self._check_file()
        if self._check_file():
            logger.info("  - GOAF filepath exists.")
            with open(self._filepath, "r") as read_content:
                temp_content = read_content.readlines()
                self._readlines = []
                for line in temp_content:
                    if not line.startswith("!") and not line.strip() == "":
                        line = line.strip()
                        line_category = self._get_go_category_from_line(line)
                        if line_category in go_categories:
                            self._readlines.append(line)
        self.terms_dict = None
        self.products_dict = None
        logger.info(f"  - GOAF created with {len(self._readlines)} annotations.")

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

    def get_all_products_for_goterm(self, goterm_id: str) -> List[str]:
        """
        This method returns all unique products associated with the GO term id.
        The return of this function is influenced by the go_categories supplied to the constructor of the GOAF!

        Args:
            goterm_id (str): a GO Term identifier, eg. GO:0003723

        Returns:
            List[str]: a List of all product (gene/gene products) gene names, eg. ['NUDT4B', ...]

        Example: for 'GO:0003723' it returns ['KMT2C', 'CLNS1A', 'ZCCHC9', 'TERT', 'BUD23', 'DDX19B', 'CCAR2', 'NAP1L4', 'SAMSN1', 'ERVK-9', 'COA6', 'RTF1', 'AHCYL1', 'SMARCA4', ... (total len = 1378)]
        """

        if self.terms_dict is None:
            self.populate_terms_dict()

        return self.terms_dict.get(goterm_id, [])

    def populate_poducts_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a product gene name and an associated GO Term.

        The result is a dictionary (self.products_dict), mapping keys (product gene names, eg. NUDT4B) to a List of all associated
        GO Terms (eg. ['GO:0003723', ...])
        """
        self.products_dict = {}
        for (
            line
        ) in (
            self._readlines
        ):  # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split("\t")
            self.products_dict.setdefault(chunks[2], set()).add(
                chunks[4]
            )  # create a key with the line's product gene name (if the key already exists, don't re-create the key - specified by the setdefault method) and add the associated GO Term to the value set. eg. {'NUDT4B': {'GO:0003723'}}, after first line is processed, {'NUDT4B': {'GO:0003723'}, 'NUDT4B': {'GO:0046872'}} after second line ...
        for (
            key,
            values,
        ) in (
            self.products_dict.items()
        ):  # the set() above prevents the value elements (GO Terms) in dictionary to be repeated
            self.products_dict[key] = list(
                values
            )  # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}

    def populate_terms_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a GO Term and it's associated product gene name.

        The result is a dictionary (self.terms_dict), mapping keys (GO Terms, eg. GO:0003723) to a List of all
        associated product gene names (eg. ['NUDT4B', ...])
        """
        self.terms_dict = {}
        for (
            line
        ) in (
            self._readlines
        ):  # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split("\t")
            self.terms_dict.setdefault(chunks[4], set()).add(
                chunks[2]
            )  # create a key with the line's GO Term (if the key already exists, don't re-create the key - specified by the setdefault method) and add the product' gene name to the value set. eg. {'GO:0003723': {'NUDT4B'}}, after first line is processed, {'GO:0003723': {'NUDT4B'}, 'GO:0046872': {'NUDT4B'}} after second line ...
        for (
            key,
            values,
        ) in (
            self.terms_dict.items()
        ):  # the previous set() prevents the value elements (product gene names) in dictionary to be repeated
            self.terms_dict[key] = list(
                values
            )  # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}

    def get_all_terms_for_product(self, product: str) -> List[str]:
        """
        Gets all GO Terms associated to a product gene name.
        The return of this function is influenced by the go_categories supplied to the constructor of the GOAF!

        Args:
          - (str) product: must be a gene name corresponding to a specific gene/gene product, eg. NUDT4B

        Returns:
          - List[str]: a List of all GO Term ids associated with the input product's gene name

        Example: for 'NUDT4B', it returns ['GO:1901911', 'GO:0071543', 'GO:0005737', 'GO:0000298', 'GO:0005634', 'GO:0034431', 'GO:0034432', 'GO:0046872', 'GO:0008486', 'GO:1901909', 'GO:0003723', 'GO:1901907', 'GO:0005829']
        """
        if self.products_dict is None:
            self.populate_poducts_dict()

        return self.products_dict.get(product, [])

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

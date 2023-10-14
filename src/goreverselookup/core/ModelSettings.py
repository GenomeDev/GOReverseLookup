from __future__ import annotations

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class ModelSettings:
    """
    Represents user-defined settings, which can be set for the model, to change the course of data processing.

      - homosapiens_only: if only homosapiens products should be queried from uniprot and ensembl # TODO: currently, this is hardcoded into requests. change this.
      - require_product_evidence_codes: # TODO implement logic
      - fisher_test_use_online_query: If True, will query the products of GO Terms (for the num_goterms_products_general inside fisher test) via an online pathway (GOApi.get_goterms).
                                      If False, fisher test will compute num_goterms_products_general (= the number of goterms associated with a product) via an offline pathway using GOAF parsing.
      - include_indirect_annotations: If True, each GO Term relevant to the analysis will hold a list of it's parents and children from the go.obo (Gene Ontology .obo) file. Also, the parents and children of GO Terms will be taken into
                                    account when performing the fisher exact test. This is because genes are annotated directly only to specific GO Terms, but they are also INDIRECTLY connected to all of the
                                    parent GO Terms, despite not being annoted directly to the parent GO Terms. The increased amount of GO Term parents indirectly associated with a gene will influence the fisher
                                    scoring for that gene - specifically, it will increate num_goterms_product_general.
                                    If False, each GO Term relevant to the analysis won't have it's parents/children computed. During fisher analysis of genes, genes will be scored only using the GO Terms that are
                                    directly annotated to the gene and not all of the indirectly associated parent GO terms.
      - datafile_paths: a dictionary that includes information about several data files used for the analysis. It has the following format:
            FILE_TYPE: {
                'organism':
                'local_filepath':
                'download_url':
            }

          example:

            'go_obo': {
                'organism': "all",
                'local_filepath': "data_files/go.obo",
                'download_url': "https://purl.obolibrary.org/obo/go.obo",
            },
            'goa_human': {
                'organism': "homo_sapiens",
                'local_filepath': "data_files/goa_human.gaf",
                'download_url': "http://geneontology.org/gene-associations/goa_human.gaf.gz"
            },
            'goa_zfin': {
                'organism': "danio_rerio",
                'local_filepath': "data_files/zfin.gaf",
                'download_url': "http://current.geneontology.org/annotations/zfin.gaf.gz"
            },
            'zfin_human_ortho_mapping': {
                'organism': "danio_rerio",
                'local_filepath': "data_files/zfin_human_ortholog_mapping.txt",
                'download_url': "https://zfin.org/downloads/human_orthos.txt"
            },
            ...
    """

    # note: specifying ModelSettings inside the ModelSettings class is allowed because of the 'from __future__ import annotations' import.
    def __init__(self) -> ModelSettings:
        self.homosapiens_only = False
        self.require_product_evidence_codes = False
        self.fisher_test_use_online_query = False
        self.include_indirect_annotations = False  # previously: include_all_goterm_parents
        self.uniprotkb_genename_online_query = False
        self.pvalue = 0.05
        self.goterms_set = []
        self.datafile_paths = {}

    @classmethod
    def from_json(cls, json_data) -> ModelSettings:
        """
        Constructs ModelSettings from a JSON representation. This is used during (ReverseLookup).load_model, when model loading
        is performed from the saved json file.
        """
        instance = cls()  # create an instance of the class
        for attr_name in dir(instance):  # iterate through class variables (ie settings in ModelSettings)
            if not callable(getattr(instance, attr_name)) and not attr_name.startswith(
                "__"
            ):
                if attr_name in json_data:  # check if attribute exists in json data
                    setattr(
                        instance, attr_name, json_data[f"{attr_name}"]
                    )  # set the attribute
                else:
                    logger.warning(
                        f"Attribute {attr_name} doesn't exist in json_data for"
                        " ModelSettings!"
                    )
        return instance

    def to_json(self):
        """
        Constructs a JSON representation of this class. This is used during (ReverseLookup).save_model to save the ModelSettings
        """
        json_data = {}
        for attr_name, attr_value in vars(self).items():
            if not callable(attr_value) and not attr_name.startswith("__"):
                json_data[attr_name] = attr_value
        return json_data

    def set_setting(self, setting_name: str, setting_value: bool):
        if hasattr(self, setting_name):
            setattr(self, setting_name, setting_value)
        else:
            logger.warning(
                f"ModelSettings has no attribute {setting_name}! Make sure to"
                " programmatically define the attribute."
            )

    def get_setting(self, setting_name: str):
        if hasattr(self, setting_name):
            return getattr(self, setting_name)
    
    def get_datafile_path(self, datafile_name:str):
        """
        Returns the local or absolute filepath to 'datafile_name'.

        Datafile names can be one of the following:
          - go_obo
          - goa_human
          - goa_zfin
          - goa_rgd
          - goa_mgi
          - goa_xenbase
          - ortho_mapping_zfin_human
          - ortho_mapping_rgd_human
          - ortho_mapping_mgi_human
          - ortho_mapping_xenbase_human
        """
        if datafile_name in self.datafile_paths:
            return self.datafile_paths[datafile_name]['local_filepath']
        else:
            logger.warning(f"datafile_name {datafile_name} not found in (ModelSettings).datafile_paths")
    
    def get_datafile_paths(self, *datafile_types:str):
        """
        Gets the paths to multiple datafiles. Datafile paths are loaded from the 'filepaths' section in input.txt.
        Allowed datafile types (the values of datafile_type parameter) are:
          - go_obo
          - goa_human
          - goa_zfin
          - goa_rgd
          - goa_mgi
          - goa_xenbase
          - ortho_mapping_zfin_human
          - ortho_mapping_rgd_human
          - ortho_mapping_mgi_human
          - ortho_mapping_xenbase_human

        Returns a dictionary with the keys corresponding to input datafile_types and the values to the found filepaths.

        Example call:
        (ReverseLookup).get_datafile_paths('go_obo', 'goa_human')
        -> returns:
            {
                'go_obo': "PATH_TO_GO_OBO_FILE" or None,
                'goa_human': "PATH_TO_GOAF_FILE" or None
            }

        Note: You can call this function with only one string parameter "ALL" in order to receive a dictionary of all possible
        datafile types.
        """
        all_keys = [
            "go_obo",
            "goa_human",
            "goa_zfin",
            "goa_rgd",
            "goa_mgi",
            "goa_xenbase",
            "ortho_mapping_zfin_human",
            "ortho_mapping_rgd_human",
            "ortho_mapping_mgi_human",
            "ortho_mapping_xenbase_human"
        ]
        if len(datafile_types) == 1 and datafile_types[0] == "ALL":
            input_keys = all_keys
        else:
            input_keys = datafile_types
        
        if any(input_key not in all_keys for input_key in input_keys):
            raise Exception(f"One or more input keys from {input_keys} are not valid. Valid keys are: {all_keys}")

        result = {}
        for datafile_type in input_keys:
            result[datafile_type] = self.get_datafile_path(datafile_type)
        return result

    def get_datafile_url(self, datafile_name:str):
        """
        Returns the download url to 'datafile_name'.

        Datafile names can be one of the following:
          - go_obo
          - goa_human
          - goa_zfin
          - goa_rgd
          - goa_mgi
          - goa_xenbase
          - ortho_mapping_zfin_human
          - ortho_mapping_rgd_human
          - ortho_mapping_mgi_human
          - ortho_mapping_xenbase_human
        """
        if datafile_name in self.datafile_paths:
            return self.datafile_paths[datafile_name]['download_url']
        else:
            logger.warning(f"datafile_name {datafile_name} not found in (ModelSettings).datafile_paths")
    
    def get_datafile_urls(self, *datafile_types:str):
        """
        Gets the download urls to multiple datafiles. Datafile urls are loaded from the 'filepaths' section in input.txt.
        Allowed datafile types (the values of datafile_types parameter) are:
          - go_obo
          - goa_human
          - goa_zfin
          - goa_rgd
          - goa_mgi
          - goa_xenbase
          - ortho_mapping_zfin_human
          - ortho_mapping_rgd_human
          - ortho_mapping_mgi_human
          - ortho_mapping_xenbase_human

        Returns a dictionary with the keys corresponding to input datafile_types and the values to the found filepaths.

        Example call:
        (ReverseLookup).get_datafile_paths('go_obo', 'goa_human')
        -> returns:
            {
                'go_obo': "GO_OBO_FILE_DOWNLOAD_URL" or None,
                'goa_human': "GOAF_FILE_DOWNLOAD_URL" or None
            }

        Note: You can call this function with only one string parameter "ALL" in order to receive a dictionary of all possible
        datafile types.
        """
        all_keys = [
            "go_obo",
            "goa_human",
            "goa_zfin",
            "goa_rgd",
            "goa_mgi",
            "goa_xenbase",
            "ortho_mapping_zfin_human",
            "ortho_mapping_rgd_human",
            "ortho_mapping_mgi_human",
            "ortho_mapping_xenbase_human"
        ]
        if len(datafile_types) == 1 and datafile_types[0] == "ALL":
            input_keys = all_keys
        else:
            input_keys = datafile_types
        
        if any(input_key not in all_keys for input_key in input_keys):
            raise Exception(f"One or more input keys from {input_keys} are not valid. Valid keys are: {all_keys}")

        result = {}
        for datafile_type in input_keys:
            result[datafile_type] = self.get_datafile_url(datafile_type)
        return result
    
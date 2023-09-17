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
    """

    # note: specifying ModelSettings inside the ModelSettings class is allowed because of the 'from __future__ import annotations' import.
    def __init__(self) -> ModelSettings:
        self.homosapiens_only = False
        self.require_product_evidence_codes = False
        self.fisher_test_use_online_query = False
        self.include_indirect_annotations = False  # previously: include_all_goterm_parents
        
        self.uniprotkb_genename_online_query = False
        self.pvalue = 0.05

    @classmethod
    def from_json(cls, json_data) -> ModelSettings:
        """
        Constructs ModelSettings from a JSON representation. This is used during (ReverseLookup).load_model, when model loading
        is performed from the saved json file.
        """
        instance = cls()  # create an instance of the class
        for attr_name in dir(
            instance
        ):  # iterate through class variables (ie settings in ModelSettings)
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

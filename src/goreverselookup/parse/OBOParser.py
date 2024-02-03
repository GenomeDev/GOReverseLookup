import networkx as nx
import os
from typing import Set, List, Dict, Optional

from ..util.Timer import Timer
from ..util.FileUtil import FileUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger

# workaround class, as importing GOTerm would cause a circular import
class GOTerm_placeholder:
    def __init__(self, id: str, name: Optional[str] = None, description: Optional[str] = None, category: Optional[str] = None, parent_term_ids: Optional[List[str]] = None, is_obsolete:bool = False, weight: float = 1.0, products: List[str] = [], http_error_codes:dict={}):
        """
        A placeholder class to construct GO Terms

        Args:
            id (str): The ID of the GO term.
            SOIs (List): The associated states of interest
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            products (list): Products associated with the term (optional).
            category (str): biological_process, molecular_activity or cellular_component
            parent_term_ids (list[str]): GO ids of the parent terms (parsed from .obo)
            is_obsolete (bool): if the term is labelled as obsolete in the .obo file
        """
        self.id = id
        self.name = name
        self.description = description
        self.products = products
        self.category = category
        self.parent_term_ids = parent_term_ids
        self.is_obsolete = is_obsolete


class OboParser:
    def __init__(self, obo_filepath: str = "data_files/go.obo", obo_download_url:str = "https://purl.obolibrary.org/obo/go.obo"):
        """
        Parses the Gene Ontology OBO file.

        Params:
          - (str) obo_filepath: the filepath to the obo file
        """
        self.format_version = None
        self.data_version = None
        self.ontology = None

        FileUtil.download_file(filepath=obo_filepath, download_url=obo_download_url)
        dag = (
            nx.MultiDiGraph()
        )  # DiGraph cannot store multiple edges, therefore use MiltiDiGraph!!!

        def _reset_term_data():
            """
            Used during file parsing
            """
            return {
                "id": None,  # obo: id
                "name": None,  # obo: name
                "category": None,  # obo: namespace
                "description": None,  # obo: definition
                "parent_term_ids": [],  # obo: is_a
                "is_obsolete": False,
            }

        # read all GO terms from the OBO file
        all_goterms = {}  # mapping of all go ids to GOTerm objects
        all_valid_goterms = {} # mapping of all valid GOTerms (which are not obsolete)
        all_valid_BP_goterms = {} # mapping of all valid "Biological Process" GO Terms
        all_valid_MF_goterms = {} # mapping of all valid "Molecular Function" GO Terms
        all_valid_CC_goterms = {} # mapping of all valid "Cellular Component" GO Terms
        all_obsolete_goterms = {} # mapping of all obsolete GOTerms

        # exe version bugfix:
        if not os.path.exists(obo_filepath):
            obo_filepath = "data_files/go.obo"

        with open(obo_filepath, "r") as obo_file:
            term_data = {
                "id": None,  # obo: id
                "name": None,  # obo: name
                "category": None,  # obo: namespace
                "description": None,  # obo: def
                "parent_term_ids": [],  # obo: is_a
                "is_obsolete": False,
            }
            is_obsolete = False
            for line in obo_file:  # also strips \n etc from lines
                line = line.strip()
                if line == "":
                    continue
                if "[Typedef]" in line:  # typedefs are at the end of the file, no more go term data is expected
                    break
                if "format-version" in line:
                    self.format_version = line.split(" ")[1]
                if "data-version" in line:
                    self.data_version = line.split(" ")[1]
                if "ontology" in line:
                    self.ontology = line.split(" ")[1]
                if "[Term]" in line:
                    if term_data["id"] is not None:  # term_data != _reset_term_data check is used so that if all values in JSON are none, the goterm creation block isn't executed
                        # 'is_obsolete' is not present in all GO Terms. If it isn't present, set 'is_obsolete' to false
                        if "is_obsolete" not in term_data:
                            term_data["is_obsolete"] = False

                        current_goterm = GOTerm_placeholder(
                            id=term_data["id"],
                            name=term_data["name"],
                            category=term_data["category"],
                            description=term_data["description"],
                            parent_term_ids=term_data["parent_term_ids"],
                            is_obsolete=term_data["is_obsolete"],
                        )
                        all_goterms[current_goterm.id] = current_goterm
                        if is_obsolete == True:
                            all_obsolete_goterms[current_goterm.id] = current_goterm
                        else:
                            all_valid_goterms[current_goterm.id] = current_goterm
                        
                        match current_goterm.category:
                            case "molecular_function":
                                all_valid_MF_goterms[current_goterm.id] = current_goterm
                            case "biological_process":
                                all_valid_BP_goterms[current_goterm.id] = current_goterm
                            case "cellular_component":
                                all_valid_CC_goterms[current_goterm.id] = current_goterm

                    term_data = _reset_term_data()  # reset term data for a new goterm
                else:  # Term is not in line -> line is GO Term data -> process term data in this block
                    chunks = line.split(": ", 1)  # split only first ": " element
                    line_identifier = chunks[0]
                    line_value = chunks[1]
                    match line_identifier:
                        case "id":
                            term_data["id"] = line_value
                        case "name":
                            term_data["name"] = line_value
                        case "def":
                            line_value = line_value.strip('"')  # definition line value contains double quotes in obo, strip them
                            term_data["description"] = line_value
                        case "namespace":
                            term_data["category"] = line_value
                        case "is_a":
                            line_value = line_value.split(" ")[0]  # GO:0000090 ! mitotic anaphase -> split into GO:0000090
                            term_data["parent_term_ids"].append(line_value)
                        case "is_obsolete":
                            is_obsolete = True if line_value == "true" else False
                            term_data["is_obsolete"] = is_obsolete

        # all go terms from OBO are now constructed as GOTerm objects in all_goterms dictionary
        # create a Direcected Acyclic Graph from the created GO Terms
        for goid, goterm in all_goterms.items():
            assert isinstance(goterm, GOTerm_placeholder)
            if dag.has_node(goterm.id) is False:
                dag.add_node(goterm.id)
            for parent_id in goterm.parent_term_ids:
                # nodes are automatically added if they are not yet in the graph when using add_edge
                dag.add_edge(all_goterms[parent_id].id, goterm.id)

        self.filepath = obo_filepath
        self.dag = dag
        self.all_goterms = all_goterms
        self.all_valid_goterms = all_valid_goterms
        self.all_valid_BP_goterms = all_valid_BP_goterms
        self.all_valid_MF_goterms = all_valid_MF_goterms
        self.all_valid_CC_goterms = all_valid_CC_goterms
        self.all_obsolete_goterms = all_obsolete_goterms
        self.previously_computed_parents_cache = {}  # cache dictionary between already computed goterms and their parents
        self.previously_computed_children_cache = {} # cache dictionary between already computed goterms and their children
        logger.info("Obo parser init completed.")
    
    def get_indirect_annotations(self, term_id:str, indirect_annotations_direction:str, max_depth:int, return_as_json:bool=False, ordered:bool=True):
        """
        Returns indirect annotations to the GO term specified by 'term_id'. 'indirect_annotations_direction' must be either "p" or "c", to specify whether
        this function returns parent or child GO terms, respectively, of the 'term_id'.
        
        Parameters:
          - (str) term_id: The GO Term whose parents you wish to obtain
          - (str) indirect_annotations_direction: Either "p" to obtain the parents or "c" to obtain the children. If None is passed (if 'p' or 'c' was not set as an optional setting to include_indirect_annotations),
                                                  then this function will default to using parents as indirect annotations
          - (int) max_depth: The maximum depth of indirect annotations. Indirect annotations beyond max depth will not be returned. If max depth isn't specified, or if max depth is set to -1, then full depth (the full tree) will be returned.
          - (bool) ordered: If True,  parents will be returned topologically (closest parents will be listed first in the returned list)
          - (bool) return_as_json
        """
        if not isinstance(max_depth, int):
            max_depth = int(max_depth)
            
        if indirect_annotations_direction == None:
            indirect_annotations_direction = "p"
            
        if indirect_annotations_direction != "p" and indirect_annotations_direction != "c":
            raise Exception(f"OboParser error: get_indirect_annotations only accepts 'p' or 'c' for indirect_annotations_direction. Passed indirect_annotations_direction: '{indirect_annotations_direction}'")
        
        indirects = None
        if indirect_annotations_direction == "p":
            indirects = self.get_parent_terms(term_id=term_id, return_as_json=return_as_json, ordered=ordered)
        if indirect_annotations_direction == "c":
            indirects = self.get_child_terms(term_id=term_id, return_as_json=return_as_json, ordered=ordered)
        
        # apply max depth
        if max_depth != -1 and max_depth is not None:
            indirects = indirects[:max_depth] # return 'up to max depth' elements
        
        return indirects

    def get_parent_terms(
        self, term_id: str, return_as_json:bool = False, ordered: bool = True
    ):
        """
        Gets all of GO Term parents of 'term_id'.

        Parameters:
          - (str) term_id: The GO Term whose parents you wish to obtain
          - (bool) return_as_json: If False, will return a list of string ids of parent GO Terms.
                                    If True, will return a list of GO Term parents represented as JSON objects, with the following structure: {'id':xxxx, 'name':xxxx, 'description':xxxxx, 'products':[...], 'category':xxxx, 'parent_term_ids':[...], 'is_obsolete':True/False}
          - (bool) ordered: If True, parents will be returned topologically (closest parents will be listed first in the returned list)

        Returns: A list of parent GO Terms (either ids or classes)
        """
        # attempt to cache old data
        if term_id in self.previously_computed_parents_cache:
            return self.previously_computed_parents_cache[term_id]

        Timer(millisecond_init=True)
        parents = []  # WARNRING!! Using set() DESTROYS THE ORDER OF PARENTS !!!
        ancestors = nx.ancestors(self.dag, term_id)
        # logger.debug(f"nx.ancestors elapsed: {timer.get_elapsed_formatted('milliseconds', reset_start_time=True)}")

        if ordered is True:
            # Calculate the distance from each ancestor to the given node
            distances = {
                ancestor: nx.shortest_path_length(
                    self.dag, source=ancestor, target=term_id
                )
                for ancestor in ancestors
            }
            # Sort ancestors by distance in ascending order
            sorted_distances = dict(sorted(distances.items(), key=lambda item: item[1]))
            ancestors = sorted_distances.keys()
        # logger.debug(f"ancestors ordering elapsed: {timer.get_elapsed_formatted('milliseconds', reset_start_time=True)}")

        for parent_id in ancestors:
            if return_as_json is True:
                parents.append(self.all_goterms[parent_id].__dict__)
            else:
                parents.append(parent_id)

        # cache
        self.previously_computed_parents_cache[term_id] = parents

        return parents

    def get_child_terms(
        self, term_id: str, return_as_json: bool = False, ordered: bool = True
    ):
        """
        Gets all of GO Term children of 'term_id'.

        Parameters:
          - (str) term_id: The GO Term whose children you wish to obtain
          - (bool) return_as_json: If False, will return a list of string ids of parent GO Terms.
                                    If True, will return a list of GO Term parents represented as JSON objects, with the following structure: {'id':xxxx, 'name':xxxx, 'description':xxxxx, 'products':[...], 'category':xxxx, 'parent_term_ids':[...], 'is_obsolete':True/False}
          - (bool) ordered: If True, parents will be returned topologically (closest children will be listed first in the returned list)

        Returns: A list of children GO Terms (either ids or classes)

        Example usage: We want to query child terms of 'GO:0003924 GTPase activity'
        When looking at Amigo2's inferred tree view for GO:0003924 (https://amigo.geneontology.org/amigo/term/GO:0003924#display-lineage-tab), we see the following structure: 
        GO:0003924 GTPase activity
          - (is_a_relation) GO:0003925 G protein activity
          - (capable_of_relation) GO:1905360 GTPase complex
          - (is_a_relation) GO:0061791 GTPase motor activity
              - (is_a_relation) GO:1990606 membrane scission GTPase motor activity
          - GO:0034260 negative regulation of GTPase activity
          - GO:0043547 positive regulation of GTPase activity
          - GO:0043087 regulation of GTPase activity
        
        We construct the following code:
            from goreverselookup import OboParser
            obo_parser = OboParser()
            child_terms = obo_parser.get_child_terms("GO:0003924")
            -> return: 3 children: ['GO:0003925', 'GO:0061791', 'GO:1990606']
        
        As we can see, the program returns only the children terms of GO:0003924 that are quatified with the
        'is_a_relation' relationship. It correctly queries the children, even the nested children.
        """
        # attempt to cache old data
        if term_id in self.previously_computed_children_cache:
            return self.previously_computed_children_cache[term_id]

        Timer(millisecond_init=True)
        children = []  # WARNRING!! Using set() DESTROYS THE ORDER OF PARENTS !!!
        descendants = nx.descendants(self.dag, term_id)
        # logger.debug(f"nx.descendants elapsed: {timer.get_elapsed_formatted('milliseconds', reset_start_time=True)}")

        if ordered is True:
            # Calculate the distance from each ancestor to the given node
            distances = {
                descendant: nx.shortest_path_length(
                    self.dag, source=term_id, target=descendant
                )
                for descendant in descendants
            }
            # Sort ancestors by distance in ascending order
            sorted_distances = dict(sorted(distances.items(), key=lambda item: item[1]))
            descendants = sorted_distances.keys()
        # logger.debug(f"descendants ordering elapsed: {timer.get_elapsed_formatted('milliseconds', reset_start_time=True)}")

        for child_id in descendants:
            if return_as_json is True:
                children.append(self.all_goterms[child_id].__dict__)
            else:
                children.append(child_id)

        # cache and return
        self.previously_computed_children_cache[term_id] = children
        return children

    def get_goterms(self, validity="valid", go_categories:list=["molecular_function", "biological_process", "cellular_component"]): # TODO: CONTINUE WITH IMPLEMENTATION OF GO CATEGORIES
        """
        Returns all GOTerms with respect to the validity scope.

        The possible 'validity' scopes are:
          - 'all': returns all GO Terms, both valid and obsolete GO Terms
          - 'valid': returns only the valid GO Terms (who don't have is_obsolete = False)
          - 'obsolete': returns only the obsolete GO Terms (who have is_obsolete = True)
        """
        if validity == "all":
            return self.all_goterms
        if validity == "obsolete":
            return self.all_obsolete_goterms
        
        if validity == "valid":
            # return valid goterms with respect to selected go_categories
            return_result = {}
            for category in go_categories:
                match category:
                    case "biological_process":
                        return_result = {**return_result, **self.all_valid_BP_goterms}
                    case "molecular_function":
                        return_result = {**return_result, **self.all_valid_MF_goterms}
                    case "cellular_component":
                        return_result = {**return_result, **self.all_valid_CC_goterms}
            return return_result
        
        logger.warning(f"get_all_goterms was called with an inappropritate validity scope '{validity}'. Possible scopes are 'all', 'valid' and 'obsolete'.")
        return None
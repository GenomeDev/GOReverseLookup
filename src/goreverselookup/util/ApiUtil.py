# Api utility classes
from .WebsiteParser import WebsiteParser
import re


class EnsemblUtil():
    def __init__(self):
        pass

    @classmethod
    def split_ensembl_id(cls, full_ensembl_id:str, feature_identifiers:list = ["E", "FM", "G", "T", "GT", "P", "R"], feature_id_web_query=True):
        """
        Ensembl ids are constructed in the following manner:

            {STABLE_ID_PREFIX}{FEATURE}{IDENTIFIER}
        
        Examples of STABLE_ID_PREFIX: "ENS" for Homo Sapiens, "ENSDAR" for danio rerio
        Ensembl features are: E (exon), FM (Ensembl protein family), G (gene), GT (gene tree), P (protein), R (regulatory feature), T (transcript)

        To obtain the features programmatically via a web scrape of https://www.ensembl.org/info/genome/stable_ids/prefixes.html,
        combine the usage of this function with WebsiteParser.ens_feature_prefixes = WebsiteParser.parse_ensembl_stable_id_feature_prefixes_table().

            feature_prefixes_id = WebsiteParser.parse_ensembl_stable_id_feature_prefixes_table() # dict between feature ids and corresponding labels
            feature_prefixes_list = [] # extract only feature ids
            for f_id, label in feature_prefixes_id.items():
                feature_prefixes_list.append(f_id)
            EnsemblUtil.split_ensembl_id(ENS_ID, feature_prefixes_list)
        
        The parameter 'feature_id_web_query' controls if the feature prefixes are queried via the above website scrape process automatically.
        This will override the 'feature_identifiers' parameter.
        
        Example: 
        EnsemblUtil.split_ensembl_id("ENSG00132")
        -> return: {'stable_id_prefix': "ENS", 'feature_prefix': "G", 'identifier': "00132"}
        """
        if isinstance(full_ensembl_id,list):
            full_ensembl_id=full_ensembl_id[0]
            
        # construct a list of feature identifiers either via a web query or via list of feature identifiers
        feature_identifiers_final = []
        if feature_id_web_query:
            feature_prefixes_dict = WebsiteParser.parse_ensembl_stable_id_feature_prefixes_table()
            for f_id, label in feature_prefixes_dict.items():
                feature_identifiers_final.append(f_id)
        else:
            feature_identifiers_final = feature_identifiers

        identifier = re.findall(r'\d+', full_ensembl_id) # use regex to extract the numbers
        identifier = ''.join(identifier) # join the extracted numbers into a single string

        ens_id_no_numbers = full_ensembl_id.rstrip(identifier)
        feature = ""
        feature_found = False

        two_letter_features = [] # example: FM, GT
        single_letter_features = [] # example: G, T, R, P, E
        for f_id in feature_identifiers_final:
            if len(f_id) == 2:
                two_letter_features.append(f_id)
            elif len(f_id) == 1:
                single_letter_features.append(f_id)

        # iterate through two_letter_features first
        for two_letter_f_id in two_letter_features:
            if two_letter_f_id == ens_id_no_numbers[-2:]: # fi the current two letter feature id is the same as the last two letters of the ens id without numbers
                feature = two_letter_f_id
                feature_found = True
                break

        # if feature was not found yet, iterate through single letter features too
        if not feature_found:
            for single_letter_f_id in single_letter_features:
                if single_letter_f_id == ens_id_no_numbers[-1:]:
                    feature = single_letter_f_id
                    feature_found = True
                    break
        
        stable_id_prefix = ens_id_no_numbers.rstrip(feature)

        return {
            'stable_id_prefix': stable_id_prefix,
            'feature_prefix': feature,
            'identifier': identifier
        }
    
    @classmethod
    def taxon_to_ensembl_label(cls, taxon_num:str, major_species_only:bool = True):
        """
        Returns an Ensembl species label corresponding to input NCBI Taxon number ('taxon_num')

        Parameters:
          - (str) taxon_num: The number of the NCBITaxon:XXXX format in string notation. Also accepts the full NCBITaxon:XXXX notation.
          - (bool) major_species_only: For some species, Ensembl notes the label in "extended" format. For example, for NCBITaxon:10090 (House mouse, Mus musculus),
                                       Ensembl would by itself return the label "mus_musculus_balbcj". If major_species_only is set to True, then only the "major" species
                                       will be returned, in this case "mus_musculus".

        Example: We want species label for Orange clownfish (NCBITaxon:161767)
            EnsemblUtil.taxon_to_ensembl_label("161767")
            -> "amphiprion_percula"
        
        """
        if ":" in taxon_num:
            taxon_num = taxon_num.split(":")[1]
        
        if isinstance(taxon_num, int):
            taxon_num = str(taxon_num)
        
        lookup_table = WebsiteParser.get_ensembl_stable_id_prefixes_table()
        if taxon_num in lookup_table:
            label = lookup_table[taxon_num]['label']

            if major_species_only == True: # mus_musculus_balbcj
                split = label.split("_")
                if len(split) >= 2:
                    label = f"{split[0]}_{split[1]}" # mus_musculus

            return label
        
        return None
    
    @classmethod
    def taxon_to_ensembl_stable_id_prefix(cls, taxon:str):
        """
        Maps a NCBITaxon (e.g. NCBITaxon:9606) to an Ensembl stable id prefix (e.g. 'ENS')
        """
        if "NCBITaxon:" in taxon:
            taxon = taxon.replace("NCBITaxon:", "")
        return WebsiteParser.ensembl_species_to_ids_to_taxons[taxon].get('stable_id_prefix')
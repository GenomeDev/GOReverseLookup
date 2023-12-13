import requests
from bs4 import BeautifulSoup
from .JsonUtil import JsonUtil
import urllib3
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger

class WebsiteParser():
    @classmethod
    def init(cls):
        cls.ensembl_species_to_ids = cls.parse_ensembl_stable_id_prefixes_table()
        cls.ensembl_stable_id_feature_prefixes = cls.parse_ensembl_stable_id_feature_prefixes_table()
        cls.ensembl_species_to_ids_to_taxons_filepath = "data_files/ens_spec_ids_taxa.json"
        cls.ensembl_species_to_ids_to_taxons = None
        try:
            cls.ensembl_species_to_ids_to_taxons = JsonUtil.load_json(cls.ensembl_species_to_ids_to_taxons_filepath)
        except (FileExistsError, FileNotFoundError):
            logger.debug(f"'data_files/ens_spec_ids_taxa.json' WAS NOT FOUND during WebsiteParser init.")
        
        if cls.ensembl_species_to_ids_to_taxons is None or cls.ensembl_species_to_ids_to_taxons == {}:
            # compute all taxon ids for species
            cls.ensembl_species_to_ids_to_taxons = cls.fetch_ensembl_taxon_numbers_for_species()
    
    @classmethod
    def get_website_content(cls, url:str):
        """
        Returns website content as a BeautifulSoup (bs4) object.
        """
        retry_strategy = Retry(
            total=3,  # Maximum number of retry attempts (including the initial request)
            backoff_factor=0.3,  # Exponential backoff factor (0.3 means 3 * (2^0.3))
            status_forcelist=[500, 502, 503, 504],  # HTTP status codes to retry on
        )

        # Create a session with the retry strategy
        session = requests.Session()
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount('http://', adapter)
        session.mount('https://', adapter)

        try:
            response = session.get(url)
            response.raise_for_status()
            html_content = response.content
            soup = BeautifulSoup(html_content, 'html.parser')
            return soup
        except requests.exceptions.RequestException as e:
            logger.warning(f"Failed to parse {url}")
            return None

    @classmethod
    def parse_ensembl_stable_id_prefixes_table(cls):
        """
        Parses the table from https://www.ensembl.org/info/genome/stable_ids/prefixes.html

        Returns a dictionary, where the keys are latin organism names (full) and the values are their
        stable ensembl id prefixes:
        
        Example:
        {
            'Equus asinus (Donkey)': "ENSEAS"
            'Vulpes vulpes (Red fox)': "ENSVVU"
            ...
        }
        """
        url = "https://www.ensembl.org/info/genome/stable_ids/prefixes.html"
        soup = cls.get_website_content(url=url)
        res = {}
        if soup:
            main_div = soup.find('div', id='static')
            tables = main_div.find_all('table')
            table = tables[1] # there are two tables -> we are searching for the second one
            # table = soup.find('table', attrs={'class':"ss autocenter", 'style':"width: 100%", 'cellpadding':"0", 'cellspacing':"0"})
            
            table_body = table.find('tbody')
            # Find all table rows ('tr') within the table
            rows = table_body.find_all('tr')
            for row in rows:
                # find all 'td'
                tds = row.find_all('td')
                if len(tds) == 2:
                    ens_prefix = tds[0].text
                    species_name = tds[1].text
                    res[species_name] = ens_prefix

        cls.ensembl_species_to_ids = res
        return res

    @classmethod
    def get_ensembl_stable_id_prefixes_table(cls):
        """
        Gets the dictionary between taxon numbers and respective Ensembl species, stable id prefixes and the species labels. Call this function to prevent constant
        re-scraping of the ensembl website! Underlying url: https://rest.ensembl.org/info/species?content-type=application/json

        The structure of the returned dictionary is:
        {
            '8083': {'name': 'platyfish', 'label': 'xiphophorus_maculatus', 'stable_id_prefix': 'ENSXMA'}, 
            '61622': {'name': 'golden snub-nosed monkey', 'label': 'rhinopithecus_roxellana', 'stable_id_prefix': 'ENSRRO'}, 
            '161767': {'name': 'orange clownfish', 'label': 'amphiprion_percula', 'stable_id_prefix': 'ENSAPE'}, 
            '31033': {'name': 'fugu', 'label': 'takifugu_rubripes', 'stable_id_prefix': 'ENSTRU'}, 
            ...
        }
        """
        if cls.ensembl_species_to_ids_to_taxons is None and cls.ensembl_species_to_ids_to_taxons != {}:
            cls.fetch_ensembl_taxon_numbers_for_species()
        return cls.ensembl_species_to_ids_to_taxons

    @classmethod
    def fetch_ensembl_taxon_numbers_for_species(cls):
        if cls.ensembl_species_to_ids is None:
            cls.parse_ensembl_stable_id_prefixes_table()

        server = "https://rest.ensembl.org"
        ext = "/info/species?"
 
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        response_json = r.json()
        species_list = response_json.get('species')
        if species_list is None:
            raise Exception(f"Failed to fetch species list!")
        
        # convert species list to dict, where keys are display names
        species_dict = {}
        for species in species_list:
            assert isinstance(species, dict)
            if species.get('display_name') is not None and species.get('taxon_id') is not None and species.get('name') is not None:
                species_dict[species.get('display_name').lower()] = species

        result_dict = {} # key = taxon id; value = {'name' = NAME, 'stable_id_prefix' = STABLE_ID_PREFIX}
        for name, stable_id_prefix in cls.ensembl_species_to_ids.items():
            # get display name
            if "(" not in name and ")" not in name:
                logger.warning(f"Species name {name} doesn't include a display name (delimited by '(' and ')') and won't be processed!")
                continue
            display_name = name.split("(")[1] # Camarhynchus parvulus (Small tree finch) -> Small tree finch)
            display_name = display_name.split(")")[0] # Small tree finch) -> Small tree finch
            display_name = display_name.lower()
            if display_name in species_dict:
                species_taxon = species_dict[display_name].get('taxon_id')
                species_label = species_dict[display_name].get('name')
                result_dict[species_taxon] = {
                    'name': display_name,
                    'label': species_label,
                    'stable_id_prefix': stable_id_prefix
                }
        
        # write to file
        if result_dict != {}:
            JsonUtil.save_json(result_dict, cls.ensembl_species_to_ids_to_taxons_filepath)
            cls.ensembl_species_to_ids_to_taxons = result_dict

        return result_dict
    
    @classmethod
    def parse_ensembl_stable_id_feature_prefixes_table(cls):
        """
        Finds the Ensembl FEATURE prefixes. On 10.25.2023, the possible feature prefixes were: E (exon), FM (Ensembl protein family), G (gene), GT (gene tree), P (protein), R (regulatory feature), T (transcript)
        Ensembl feature prefixes are parsed from https://www.ensembl.org/info/genome/stable_ids/prefixes.html
        
        Return a dictionary, with feature prefixes as keys and labels as values:
        {
            'E': "exon",
            'FM': "Ensembl protein family"
            ...
        }
        """
        if hasattr(cls, "ensembl_stable_id_feature_prefixes"):
            if cls.ensembl_stable_id_feature_prefixes is not None and cls.ensembl_stable_id_feature_prefixes != {}:
                return cls.ensembl_stable_id_feature_prefixes
        
        url = "https://www.ensembl.org/info/genome/stable_ids/prefixes.html"
        soup = cls.get_website_content(url=url)
        res = {}
        if soup:
            main_div = soup.find('div', id='static')
            tables = main_div.find_all('table')
            table = tables[0] # there are two tables -> we are searching for the first one
            # table_body = table.find('tbody')
            # rows = table_body.find_all('tr') # Find all table rows ('tr') within the table
            rows = table.find_all('tr')
            for row in rows:
                th = row.find('th')
                td = row.find('td') # find all 'td'
                ens_prefix = th.text
                label = td.text
                res[ens_prefix] = label

        cls.ensembl_stable_id_feature_prefixes = res
        return res

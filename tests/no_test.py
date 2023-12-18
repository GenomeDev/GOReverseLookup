import os
from goreverselookup import Cacher
from goreverselookup import ReverseLookup
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test
from goreverselookup import LogConfigLoader
from goreverselookup import GOTerm
from goreverselookup import OboParser
from goreverselookup import JsonUtil

# setup logger
import logging

# logger = logging.getLogger(__name__)
# logger.addHandler(logging.StreamHandler)
LogConfigLoader.setup_logging_config(log_config_json_filepath="logging_config.json")
logger = logging.getLogger(__name__)

root_term = GOTerm(id="GO:0003924")
obo_parser = OboParser()

#root_term_children = obo_parser.get_child_terms(root_term.id)
#logger.info(f"{len(root_term_children)} children: {root_term_children}")

# print statistically relevant genes

stat_relev_genes_direct = JsonUtil.load_json("results/run7 - without indirect annotations, p=0.05, 39 genes/statistically_relevant_genes.json")
gene_names_direct = []
for gene in stat_relev_genes_direct["chronic_inflammation+:cancer+"]:
	gene_names_direct.append(gene["genename"])
print(gene_names_direct)

stat_relev_genes_indirect = JsonUtil.load_json("results/run8 - with indirect annotations, p=0.05, 130 genes/statistically_relevant_genes.json")
gene_names_indirect = []
for gene in stat_relev_genes_indirect["chronic_inflammation+:cancer+"]:
	gene_names_indirect.append(gene["genename"])
print(gene_names_indirect)

# find diff
diff = []
for gene_direct in gene_names_direct:
	if gene_direct not in gene_names_indirect:
		diff.append(gene_direct)
print(f"{len(diff)} 'direct' genes weren't found among 'indirect' genes. List: {diff}")



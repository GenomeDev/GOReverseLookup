import requests

class GProfiler:
	def __init__(self):
		pass

	@classmethod
	def convert_ncbitaxon_to_gprofiler(cls, taxon:str):
		"""
		Converts an NCBI-type taxon to a respective GProfiler taxon.
		Note: gprofiler - https://biit.cs.ut.ee/gprofiler/
		"""
		r = requests.get("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list")
		taxon_equivalents = {}
		data = r.json()
		for d in data:
			taxon_equivalents[d['taxonomy_id']] = d['id']
		return taxon_equivalents.get(taxon, None)
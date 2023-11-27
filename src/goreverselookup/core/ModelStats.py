# Stats class offers insights into the results of various model querying results.

class ModelStats:
    @classmethod
    def init(cls):
        cls.goterm_product_query_results = {}
        cls.goterm_count = 0
        cls.product_count = 0
        cls.product_ortholog_query_results = {}
        cls.product_idmapping_results = {}
        cls.product_ortholog_gOrth_query_results = {}
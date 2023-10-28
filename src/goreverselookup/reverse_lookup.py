from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from .associations import Annotations
    from .ontology import GODag

from .multiple_testing import multiple_correction
from .pvalcalc import pvalue_calculate


class ReverseLookupRecord(object):
    """Represents one result (from a single product) in the ReverseLookupStudy"""

    def __init__(
        self,
        objid,
        name=None,
        pvals={},
        study_items=set(),
        population_items=set(),
        ratio_in_study=(0, 0),
        ratio_in_pop=(0, 0),
    ):
        self.object_id = objid
        self.name = name
        self.pvals = pvals
        self.study_items = study_items
        self.population_items = population_items
        # Ex: ratio_in_pop ratio_in_study study_items p_uncorrected pop_items
        self.study_count = ratio_in_study[0]
        self.study_n = ratio_in_study[1]
        self.pop_count = ratio_in_pop[0]
        self.pop_n = ratio_in_pop[1]

    def add_pval(self, method, pvalue: float):
        """Pvalue to dict."""
        self.pvals[method] = pvalue


class GOReverseLookupStudy:
    """Runs pvalue test, as well as multiple corrections"""

    def __init__(
        self,
        anno: Annotations,  # this is annotation object. This is the population. (preprocess it to add orthologs or to propagate associations to parents). NOTE: species you add to the association object affect the result; only include the target species and the ones ortologs were searched for.
        obo_dag: GODag,  # check if it is needed?
        alpha=0.05,
        pvalcalc="fisher_scipy_stats",
        methods=None,
    ):
        self.anno = anno
        self.obo_dag = obo_dag
        self.alpha = alpha
        self.methods = methods
        if methods is None:
            self.methods = ["bonferroni"]  # add statsmodel multipletest
        self.pval_method = pvalcalc

    def run_study(
        self, studyset: Union[set[str], list[str]], **kws
    ) -> list[ReverseLookupRecord]:
        """_summary_

        Args:
            study (_type_): list of all goterms (term_id) for a process-

        Returns:
            List[ReverseLookupRecord]: _description_
        """
        """Run Gene Ontology Reverse Lookup Study"""

        if len(studyset) == 0:
            return []

        # process kwargs
        methods = kws.get("methods", self.methods)
        alpha = kws.get("alpha", self.alpha)

        # calculate the uncorrected pvalues using the pvalcalc of choice
        results = self.get_pval_uncorr(
            studyset
        )  # results is a list of ReverseLookupRecord objects
        if not results:
            return []

        # do multipletest corrections on uncorrected pvalues, add to ReverseLookupRecord objects
        self._run_multitest_corr(results, methods, alpha)

        # 'keep_if' can be used to keep only significant GO terms. Example:
        #     >>> keep_if = lambda nt: nt.p_fdr_bh < 0.05 # if results are significant
        #     >>> goea_results = goeaobj.run_study(geneids_study, keep_if=keep_if)
        if "keep_if" in kws:
            keep_if = kws["keep_if"]
            results = [r for r in results if keep_if(r)]

        # Default sort order:
        # results.sort(key=lambda r: [r.pvals["uncorrected"]])

        return results  # list of ReverseLookupRecord objects

    def get_pval_uncorr(
        self, studyset: Union[set[str], list[str]]
    ) -> list[ReverseLookupRecord]:
        """Calculate the uncorrected pvalues for study items."""
        results = []

        dict_by_object_id = self.anno.dict_from_attr("object_id")
        dict_by_term_id = self.anno.dict_from_attr("term_id")

        study2annoobjid = (
            set()
        )  # list of all annotation objects id from goterms in study
        for term_id in studyset:
            for annoobj in dict_by_term_id.get(term_id, set()):
                study2annoobjid.add(annoobj.object_id)

        for object_id in study2annoobjid:
            # for each object id (product id) calculate pvalue
            study_items = set(
                anno_obj.term_id
                for anno_obj in dict_by_object_id[object_id]
                if anno_obj.term_id in studyset
            )
            study_count = len(
                study_items
            )  # for each object id (product id) check how many goterms in study are associated to it

            if study_count == 0:
                pass  # this ensures pvalcalc is only done with annotations with at least one association to study set of goterms.

            study_n = len(studyset)  # N of study set

            population_items = set(
                anno_obj.term_id for anno_obj in dict_by_object_id[object_id]
            )
            pop_count = len(
                population_items
            )  # total number of goterms an objectid (product id) is associated in the whole population set

            pop_n = len(self.obo_dag)  # total number of goterms in population set

            one_record = ReverseLookupRecord(
                object_id,
                pvals={
                    "uncorrected": pvalue_calculate(
                        study_count, study_n, pop_count, pop_n, self.pval_method
                    )
                },
                study_items=study_items,
                population_items=population_items,
                ratio_in_study=(study_count, study_n),
                ratio_in_pop=(pop_count, pop_n),
            )

            results.append(one_record)

        return results

    def _run_multitest_corr(
        self, results: list[ReverseLookupRecord], methods: str, a: float
    ):
        pvals = [r.pvals["uncorrected"] for r in results]
        for method in methods:
            corrected_pvals = multiple_correction(pvals, method, a)
            self._update_pvalcorr(results, method, corrected_pvals)

    @staticmethod
    def _update_pvalcorr(
        results: list[ReverseLookupRecord], method: str, corrected_pvals: list[float]
    ):
        """Add data members to store multiple test corrections."""
        if corrected_pvals is None:
            return
        for rec, val in zip(results, corrected_pvals):
            rec.add_pval(method, val)


def results_intersection(
    *lists: list[ReverseLookupRecord],
) -> dict[str, list[ReverseLookupRecord]]:
    intersection_dict = defaultdict(list)
    # Create a dictionary of object_ids and their occurrences
    object_id_counts: dict[str, int] = defaultdict(int)
    for lst in lists:
        for record in lst:
            object_id_counts[record.object_id] += 1

    # Find object_ids that occur in all lists
    num_lists = len(lists)
    for object_id, count in object_id_counts.items():
        if count == num_lists:
            for lst in lists:
                for record in lst:
                    if record.object_id == object_id:
                        intersection_dict[object_id].append(record)

    return dict(intersection_dict)

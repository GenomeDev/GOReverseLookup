from __future__ import annotations
from typing import TYPE_CHECKING, Dict
from scipy.stats import binomtest, fisher_exact
import os
import gzip
import requests
import time
import json
from typing import List

if TYPE_CHECKING:
    from ..Model import ReverseLookup, Product, miRNA
    from ..parse.GOAFParser import GOAnnotationsFile

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class Metrics:
    """
    A super-class (interface) for the scoring used in the GO Reverse Lookup. It has to be implemented by a subclass, specifically you have
    to implement the 'metric' function. It is designed to measure the regulatory effect of Gene Ontology (GO) terms on a product.
    """

    def __init__(self, model: ReverseLookup):
        self.reverse_lookup = model
        self.name: str = None

    def metric(self, product: Product | miRNA):
        """
        The 'metric' function should be implemented in the subclasses of this interface.
        """
        raise NotImplementedError("Subclasses must implement metric()")


class adv_product_score(Metrics):
    """
    An advanced scoring algorithm, an implementation of the Metrics interface. It takes in  a model (in the form of a ReverseLookup object)
    and several parameters (a, b1, b2, c1, c2) in it's constructor, which are used to tune the weights given to different aspects of the scoring algorithm.

    Parameters:
      - (ReverseLookup) model: an instance of the ReverseLookup model.
      - (float) a: is used to give a base score to a product when all target processes are regulated in the same direction as the GOTerms in the list.
      - (float) b1, b2: are used to calculate the score based on the number of processes in target_processes that are regulated by GOTerms in the same (b1) or opposite (b2) direction as defined in the list.
      - (float) c1, c2: are used to adjust the score based on the number of GOTerms with direction "0"

    Scoring algorithm (explanation of the metric function):
      1. Start with score = 0.0

      2. (a) If all GO Terms of a Product instance regulate the processes of the ReverseLookup instance (eg. angio, diabetes) positively (and none negatively), then add 'a' to score. The 'direction' value of a positive regulatory term is '+', whereas direction for a negative regulatory term is '-'.

      3. (b1) For each of the processes, compute sum(goterm.weight) ** b2, for every GO Term of the product, which positively regulates the process.
         Final equation ADDED to the score is sum(b1 * sum(goterm.weight) ** b2). The first 'sum' is the sum of processes, whereas the second 'sum' is the sum of GO Terms, which pass the positive regulation check for the current process.

      4. (b2) For each of the process, compute sum(goterm.weight ** b2), for every GO Term of the product, which negatively regulates the process.
         Final equation SUBTRACTED from the score is sum(b1 * sum(goterm.weight ** b2)). The first 'sum' is the sum of processes, the second 'sum' is the sum of GO Terms, which pass the negative regulation check for the current process.

      5. (c1, c2): For the GO Terms of the product with "general" regulation (direction = 0), add a base score of c1 and add sum(c2 * goterm.weight) for every GO Term with direction = 0 (general regulation).
         Final equation ADDED to the score is score * (c1 + sum(c2 * goterm.weight)), the sum relating to the number of GO Terms with direction == 0.

    Example of calling and usage:
    1. Construct a ReverseLookup model
    model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

    2. Create the scoring object
    adv_score = adv_product_score(model)

    3. Score the products
    model.score_products(adv_score)
    """

    def __init__(
        self,
        model: ReverseLookup,
        a: float = 10,
        b1: float = 2,
        b2: float = 0.5,
        c1: float = 1,
        c2: float = 0.1,
    ):
        super().__init__(model)
        self.name = "adv_score"
        self.a = a
        self.b1 = b1
        self.b2 = b2
        self.c1 = c1
        self.c2 = c2

    def metric(self, product: Product) -> float:
        """
        An implementation of the scoring algorithm for an input Product instance.

        Parameters:
          - (Product) product: an instance of a Product

        Returns:
          - (float) score: a score according to this class' scoring algorithm.
        """
        # a list of GO Terms associated with the current Product
        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        score = 0.0

        def _opposite_direction(direction: str) -> str:
            if direction == "0":
                return "0"
            elif direction == "+":
                return "-"
            elif direction == "-":
                return "+"

        """ These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
        # Check if all target processes are regulated in the same direction as the GOTerms in the list
        # and none of them are regulated in the opposite direction
        if (
            # Check if all processes in target_processes of the ReverseLookup model
            # have a GOTerm in goterms_list that regulates it (them) in the same direction
            all(
                any(any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            and not any(
                any(any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            
        ):
            # If all target processes are regulated in the same direction, add a points to the score
            score += self.a
        """

        # Check if all target processes are regulated in the same direction as the GOTerms in the list
        # and none of them are regulated in the opposite direction
        if (
            # Check if all processes in target_processes of the ReverseLookup model
            # have a GOTerm in goterms_list that regulates it (them) in the same direction
            all(
                any(
                    any(
                        process["direction"] == goterm_process["direction"]
                        and process["process"] == goterm_process["process"]
                        for goterm_process in goterm.processes
                    )
                    for goterm in goterms_list
                )
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            and not any(
                any(
                    any(
                        _opposite_direction(process["direction"])
                        == goterm_process["direction"]
                        and process["process"] == goterm_process["process"]
                        for goterm_process in goterm.processes
                    )
                    for goterm in goterms_list
                )
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the same direction, add a points to the score
            score += self.a

        """ These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
        # Check if all target processes are regulated in the opposite direction as the GOTerms in the list
        # and none of them are regulated in the same direction
        if (
            # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            all(
                any(any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            and not any(
                any(any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the opposite direction, subtract a points from the score
            score -= self.a
        """
        # Check if all target processes are regulated in the opposite direction as the GOTerms in the list
        # and none of them are regulated in the same direction
        if (
            # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            all(
                any(
                    any(
                        _opposite_direction(process["direction"])
                        == goterm_process["direction"]
                        and process["process"] == goterm_process["process"]
                        for goterm_process in goterm.processes
                    )
                    for goterm in goterms_list
                )
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            and not any(
                any(
                    any(
                        process["direction"] == goterm_process["direction"]
                        and process["process"] == goterm_process["process"]
                        for goterm_process in goterm.processes
                    )
                    for goterm in goterms_list
                )
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the opposite direction, subtract a points from the score
            score -= self.a

        """ These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the same direction as defined in the list
        score += sum(
            (self.b1 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes)
                ) ** self.b2)
                for process in self.reverse_lookup.target_processes
        )
        """

        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the same direction as defined in the list
        sum_weights = 0
        for goterm in goterms_list:
            for goterm_process in goterm.processes:
                for process in self.reverse_lookup.target_processes:
                    if (
                        goterm_process["direction"] == process["direction"]
                        and goterm_process["process"] == process["process"]
                    ):
                        sum_weights += goterm.weight

                score += self.b1 * sum_weights**self.b2
                sum_weights = 0

        """ These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the oposite direction as defined in the list
        score -= sum(
            (self.b1 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes)
                ) ** self.b2)
                for process in self.reverse_lookup.target_processes
        )
        """

        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the oposite direction as defined in the list
        sum_weights = 0
        for goterm in goterms_list:
            for goterm_process in goterm.processes:
                for process in self.reverse_lookup.target_processes:
                    if (
                        goterm_process["direction"]
                        == _opposite_direction(process["direction"])
                        and goterm_process["process"] == process["process"]
                    ):
                        sum_weights += goterm.weight

                score -= self.b1 * sum_weights**self.b2
                sum_weights = 0

        """
        # These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
        # Calculate the score by multiplying the current score with a factor based on the number of GOTerms with direction "0"
        score = score * (
            self.c1  # Start with a base factor of 1
            + (self.c2  # Add a factor based on a constant value c
                * sum(  # Multiply c by the sum of weights of all GOTerms with direction "0"
                    goterm.weight  # Get the weight of each GOTerm
                    for goterm in goterms_list  # Iterate over all GOTerms in the list
                    if any(p['direction'] == 0 for p in goterm.processes)  # Only consider GOTerms with direction "0"
                )
            )
        )
        """

        # Calculate the score by multiplying the current score with a factor based on the number of GOTerms with direction "0"
        sum_weights = 0
        for goterm in goterms_list:
            for goterm_process in goterm.processes:
                if goterm_process["direction"] == "0":
                    sum_weights += goterm.weight
        score *= self.c1 + self.c2 * sum_weights
        sum_weights = 0

        return score


class nterms(Metrics):
    """
    An implementation of the Metrics interface, it scores the products by positive, negative or general regulation of a speciffic process.
    It takes in a model (in the form of a ReverseLookup object) in it's constructor.

    Parameters:
      - (ReverseLookup) model: an instance of the ReverseLookup model.

    Scoring algorithm (explanation of the metric function):
      Create an empty nterms_dict, where descriptive regulatory keys (eg. 'angio+', 'angio-', 'angio0') will be mapped to a count of terms regulating a specific process in a specific direction
      For each process in the model's (ReverseLookup) 'target_processes':
        a) create the following keys in nterms_dict: '{process}+', '{process}-', '{process}0'; if process is 'angio', then 'angio+', 'angio-', 'angio0' will be the keys in nterms_dict
        b) populate each of the keys with the count of GO Terms, which positively (direction == '+'), negatively (direction == '-') or generally (direction == '0') regulate the process

    Example of calling and usage:
    1. Create a ReverseLookup model
    model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

    2. Create the scoring object
    nterms_score = nterms(model)

    3. Use the scoring object using model.score_products
    model.score_products(nterms_score)
    """

    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "nterms"

    def metric(self, product: Product) -> dict:
        """
        An implementation of the scoring algorithm for an input Product instance.

        Parameters:
          - (Product) product: an instance of a Product

        Returns:
          - (dict) nterms_dict: a dictionary with (ReverseLookup).target_processes * 3 keys. Each process of a ReverseLookup instance has 3 keys,
                                '{process}+', '{process}-', '{process}0'. Each key has an integer count value of the amount of GO Terms of the input Product instance,
                                which positively (direction == '+'), negatively (direction == '-') or generally (direction == '0') regulate a speciffic process.

                                For a ReverseLookup model with defined processed 'angio' and 'diabetes', the returned dictionary would have 6 keys:
                                angio+, angio-, angio0, diabetes+, diabetes-, diabetes0
        """
        # A list of GO Terms associated with the current Product
        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        # An empty dictionary to store the count of GOTerms for each process and direction
        nterms_dict = {}

        # Iterate over each process in the target_processes list
        for process in self.reverse_lookup.target_processes:
            # Count the number of GOTerms that have a direction of "+" and a process matching the current process
            # nterms_dict[f"{process['process']}+"] = sum(1 for goterm in goterms_list if (goterm.direction == "+" and process['process'] == goterm.process)) # These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
            nterms_dict[f"{process['process']}+"] = sum(
                1
                for goterm in goterms_list
                if any(
                    goterm_process["direction"] == "+"
                    and process["process"] == goterm_process["process"]
                    for goterm_process in goterm.processes
                )
            )

            # Count the number of GOTerms that have a direction of "-" and a process matching the current process
            # nterms_dict[f"{process['process']}-"] = sum(1 for goterm in goterms_list if (goterm.direction == "-" and process['process'] == goterm.process)) # These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
            nterms_dict[f"{process['process']}-"] = sum(
                1
                for goterm in goterms_list
                if any(
                    goterm_process["direction"] == "-"
                    and process["process"] == goterm_process["process"]
                    for goterm_process in goterm.processes
                )
            )

            # Count the number of GOTerms that have a direction of "0" and a process matching the current process
            # nterms_dict[f"{process['process']}0"] = sum(1 for goterm in goterms_list if (goterm.direction == "0" and process['process'] == goterm.process)) # These scorings worked, when "direction" was a direct submember of a GOTerm instance. Now, each GOTerm has a list of dictionaries, each dictionary representing one process with keys "process" and "direction"
            nterms_dict[f"{process['process']}0"] = sum(
                1
                for goterm in goterms_list
                if any(
                    goterm_process["direction"] == "0"
                    and process["process"] == goterm_process["process"]
                    for goterm_process in goterm.processes
                )
            )

        # Return the dictionary containing the count of GOTerms for each process and direction
        return nterms_dict


class binomial_test(Metrics):
    def __init__(self, model: ReverseLookup, goaf: GOAnnotationsFile = None):
        super().__init__(model)
        if goaf is None:
            self.goaf = model.goaf
        else:
            self.goaf = goaf
        self.name = "binomial_test"
        self._num_all_goterms = 0

    def metric(self, product: Product) -> Dict:
        # get the count of all GO terms from the GOAF
        if self._num_all_goterms == 0:
            self._num_all_goterms = len(self.goaf.get_all_terms())

        results_dict = {}

        for process in self.reverse_lookup.target_processes:
            process_goterms_list = self.reverse_lookup.get_all_goterms_for_process(
                process["process"]
            )  # get all (positive, negative, neutral) GO terms for this process from the input file
            num_goterms_product_general = len(
                self.goaf.get_all_terms_for_product(product.genename)
            )  # get all GO terms associated with this product from the GOAF
            num_goterms_all_general = self._num_all_goterms
            for direction in ["+", "-"]:
                # num goterms associated with input Product p AND the current process (including process direction)
                num_goterms_product_process = sum(
                    1
                    for goterm in process_goterms_list
                    if (
                        any(i["direction"] == direction for i in goterm.processes)
                        and (
                            any(
                                product_id in goterm.products
                                for product_id in product.id_synonyms
                            )
                            or product.genename in goterm.products
                        )
                    )
                )
                # num goterms associated with this process (incl. direction)
                num_goterms_all_process = sum(
                    1
                    for goterm in process_goterms_list
                    if any(i["direction"] == direction for i in goterm.processes)
                )

                # time for Binomial test and "risk ratio"
                # binom = binomtest(num_goterms_product_process, num_goterms_all_process,
                #                  (num_goterms_product_general/num_goterms_all_general), alternative='greater')

                binom = binomtest(
                    num_goterms_product_process,
                    num_goterms_all_process,
                    (
                        num_goterms_product_general / num_goterms_all_general
                        if num_goterms_all_general != 0
                        else 0
                    ),
                    alternative="greater",
                )  # bugfix: ZeroDivisionError
                binom_pvalue = binom.pvalue

                if (
                    num_goterms_product_general != 0 and num_goterms_all_general != 0
                ):  # bugfix: ZeroDivisionError
                    risk_ratio = (
                        num_goterms_product_process / num_goterms_all_process
                        if num_goterms_all_process != 0
                        else 0
                    ) / (num_goterms_product_general / num_goterms_all_general)
                else:
                    risk_ratio = 0

                fold_enrichment_score = 0
                if (
                    num_goterms_all_process != 0
                    and num_goterms_product_general != 0
                    and num_goterms_all_general != 0
                ):
                    fold_enrichment_score = num_goterms_product_process / (
                        num_goterms_all_process
                        * (num_goterms_product_general / num_goterms_all_general)
                    )

                results_dict[f"{process['process']}{direction}"] = {
                    # "n_prod_process" : num_goterms_product_process,
                    # "n_all_process" : num_goterms_all_process,
                    # "n_prod_general" : num_goterms_product_general,
                    # "n_all_general" : num_goterms_all_general,
                    "num": num_goterms_product_process,
                    "expected": num_goterms_all_process
                    * (
                        num_goterms_product_general / num_goterms_all_general
                        if num_goterms_all_general != 0
                        else 0
                    ),
                    "fold_enrichment": (
                        fold_enrichment_score
                    ),  # bugfix: ZeroDivisionError
                    "pvalue": binom_pvalue,
                    "risk_ratio": risk_ratio,
                }

        # all_target_pvalues = [results_dict[f"{process['process']}{process['direction']}"]['pvalue'] for process in self.reverse_lookup.target_processes]
        # combined_p = combine_pvalues(all_target_pvalues)
        # results_dict["comb_binom_pvalue"] = combined_p.pvalue
        # combined_rr = statistics.mean([results_dict[f"{process['process']}{process['direction']}"]['risk_ratio'] for process in self.reverse_lookup.target_processes])
        # results_dict["comb_risk_ratio"] = combined_rr

        return results_dict


class fisher_exact_test(Metrics):
    """
    Fisher exact test.

    Parameters:
      - (ReverseLookup) model: an instance of the ReverseLookup model
      - (GOAnnotationsFile) goaf: an instance of the GOAnnotationsFile. If it isn't supplied, the GOAnnotationsFile from model will be used.
                                  It is advisable to not supply a custom GOAF file, since the GOAF is taken from the model.

    Note: If the setting reverse_lookup.model_settings.fisher_test_use_online_query is True (inside ModelSettings in ReverseLookup),
    then num_goterms_product_general will be determined via an https query. Otherwise, num_goterms_product_general will be determined from
    the GOAF.

    Each process can have 2 sets of GO terms - one set includes GO terms, which promote,
    and the other set includes GO terms which inhibit the process. The “general” set
    contains all GO terms in existence (found in the GOAF). Note that when supplying a ReverseLookup
    and a GOAnnotationsFile instances to the constructor, ideally they should match in their go_categories.
    If ReverseLookup and GOAF don't match in go_categories, then GOAF will be recalculated using ReverseLookup's
    go_categories.

    For each gene, we construct a contingency table and calculate the p values according
    to Fischer’s exact test.

    The gene is a candidate gene for a “positive” (stimulatory) cross section
    (eg angiogenesis + diabetes) only if p<0.05 for all stimulatory processes
    (angio+, diab+) and p>0.05 for all inhibitory processes (angio-,diab-).
    Because we calculate 2*n_processes (each process has + or - direction) p values
    for each gene, we need to calculate the final overall p value using BH correction.
    Note that BH correction is calculated in (ReverseLookup).score_products after all products have
    been scored.

    Example: consider process “angiogenesis+”. Let there exist 100 GO terms, which
    stimulate angiogenesis. The gene in question is SOX2, which is associated in 10
    of the 100 GO terms of angiogenesis+. Gene SOX2 is also associated with 100 GO
    terms in the “general” set (containing all existing GO terms - 10000 GO terms).

                                    | n GOt (contains SOX2) | n GOt (doesnt contain SOX2)  | total
    --------------------------------------------------------------------------------------------------
    set of GOt for angio+       | 10	                | 100-10=90	                   | 100
    --------------------------------------------------------------------------------------------------
    general set of GOt (all)    |                       |                              |
    MINUS set of GOt for angio+ | (100-10)=90	        | (9900-(100-10))=9810	       | 9900
    --------------------------------------------------------------------------------------------------
    total	                    | 100	                | 9900	                       | 10000

    The same table with filled out base code variables for the Fisher's test:
                                | n GOt (contains SOX2)         | n GOt (doesnt contain SOX2)  | total
    --------------------------------------------------------------------------------------------------
    set of GOt for angio+       | num_goterms_product_process   | ?                            | num_goterms_all_process
    --------------------------------------------------------------------------------------------------
    general set of GOt (all)    |                               |                              |
    MINUS set of GOt for angio+ | ?                             | ?                            | num_goterms_all_general
    --------------------------------------------------------------------------------------------------
    total                       | num_goterms_product_general   | num_goterms_all_general      |

    The complete table can now be calculated:
                                    | n GOt (contains SOX2)         | n GOt (doesnt contain SOX2)                               | total
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    set of GOt for angio+           | num_goterms_product_process   | num_goterms_all_process - num_goterms_product_process     |  num_goterms_all_process
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    general set of GOt (all)        | num_goterms_product_general - | num_goterms_all_general -                                 |
    MINUS set of GOt for angio+     | num_goterms_product_process   | (num_goterms_all_process - num_goterms_product_process)   |
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    total                           | num_goterms_product_general   | num_goterms_all_general                                   |


    Ladi original:
    Recimo da imaš 3 processe, katere zelis da se skupaj zgodijo: za vsak process določiš
    dva seta, enega z gotermi, ki sprožajo proces in enega z gotermi, ki zavirajo process
    (negative regulation + posredno zaviralci prek nasprotujočih se procesov).
    Nato določiš še splošni set (vsebuje vse goterme, ki se pojavljajo v tvojem tarčnem
    "organizmu", to ubistu pomeni kar vsi gotermi v GO, izjemno redko bi ožal, če bi imel
    kak specifičen primer), s katerim boš primerjal vsakega od prej definiranih dveh setov
    za vsak proces. Nato za vsak gen "nardis" kontingencno tabelco in izračunaš p-vrednost
    po Fischer's exact test in s tem dobimo p vrednost.

    In gen je kandidatni samo ce je p vrednost recmo manjsa od 0.05 v vseh procesih in ni
    manjsa od 0.05 v nobenem od setov gotermov, ki zavirajo proces. Aha pa se to:
    ker ubistu vsakemu genu izracunava 2*N(processov) p-vrednosti, je to treba upostevat,
    recimo bonferoni popravek za p vrednost je (target p-value)/number of test. No cisto
    nakoncu lahko iz Nprocessov p-vrednosti za sete, ki pospesujejo process izracunava neko
    meta-p vrednost, ki bi bla recmo povprecje al pa neki tazga (to morm se preucit), in nato
    bi za razvrstitev kandidatnih genov zracunala se -log(p-value) in risk ratio.
    """

    def __init__(self, model: ReverseLookup, goaf: GOAnnotationsFile = None):
        super().__init__(model)
        if goaf is None:
            self.goaf = model.goaf
        else:
            self.goaf = goaf
        self.name = "fisher_test"
        self._num_all_goterms = 0
        if self.reverse_lookup.model_settings.fisher_test_use_online_query is True:
            self.online_query_api = self.reverse_lookup.go_api
        if self.goaf.go_categories != self.reverse_lookup.go_categories:
            logger.warning(
                "GOAF categories don't match ReverseLookup model GO categories!"
            )
            logger.warning(f"  - GOAF GO categories: {self.goaf.go_categories}")
            logger.warning(
                f"  - ReverseLookup GO categories: {self.reverse_lookup.go_categories}"
            )
            logger.info(
                "GOAF will be recalculated using the ReverseLookup's GO categories:"
                f" {self.reverse_lookup.go_categories}"
            )
            self.goaf = GOAnnotationsFile(
                filepath=self.reverse_lookup.datafile_paths["goaf_filepath"],
                go_categories=self.reverse_lookup.go_categories,
            )

    def metric(self, product: Product) -> Dict:
        D_DEBUG_CALCULATE_DESIRED_N_PROD_PROCESS = True  # TODO: delete this # calculates num_goterms_product_process which would be sufficient for the product's statistical importance (p < 0.05)

        if self._num_all_goterms == 0:
            self._num_all_goterms = len(self.goaf.get_all_terms())

        results_dict = {}

        for process in self.reverse_lookup.target_processes:  # example self.reverse_lookup.target_processes: [0]: {'process': 'chronic_inflammation', 'direction': '+'}, [1]: {{'process': 'cancer', 'direction': '+'}}
            process_goterms_list = self.reverse_lookup.get_all_goterms_for_process(process["process"])  # all GO Term ids associated with a specific process (eg. angio, diabetes, obesity) for the current MODEL
            num_goterms_all_general = (self._num_all_goterms)  # number of all GO Terms from the GO Annotations File (current total: 18880)

            # num_goterms_product_general ... # all GO Terms associated with the current input Product instance (genename) from the GO Annotation File
            #   - can be queried either via online or offline pathway (determined by model_settings.fisher_test_use_online_query)
            #   - can have all parent terms (indirectly associated terms) added to the count (besides only directly associated GO terms) - determined by model_settings.include_indirect_annotations
            if self.reverse_lookup.model_settings.fisher_test_use_online_query is True:  # online pathway: get goterms associated with this product via a web query
                goterms_product_general = self.online_query_api.get_goterms(product.uniprot_id, go_categories=self.reverse_lookup.go_categories)
                if goterms_product_general is not None:
                    num_goterms_product_general = len(goterms_product_general)
                else:
                    # skip iteration, there was an error with querying goterms associated with a product
                    logger.warning(f"Online query for GO Terms associated with {product.uniprot_id} failed! Product: {json.dumps(product.__dict__)}")
                    continue
                # num_goterms_product_general = len(self.online_query_api.get_goterms(product.uniprot_id, go_categories=self.reverse_lookup.go_categories))
                logger.debug(f"Fisher test online num_goterms_product_general query: {num_goterms_product_general}")
            else:  # offline pathway: get goterms from GOAF
                if self.reverse_lookup.model_settings.include_indirect_annotations == True:
                    goterms_product_general = self.goaf.get_all_terms_for_product(product.genename, indirect_annotations=True, obo_parser=self.reverse_lookup.obo_parser)
                else:
                    goterms_product_general = self.goaf.get_all_terms_for_product(product.genename)
                num_goterms_product_general = len(goterms_product_general)  # all GO Terms associated with the current input Product instance (genename) from the GO Annotation File

            # find the number of parent indirect annotations
            if self.reverse_lookup.model_settings.include_indirect_annotations is True:
                # include all parent goterms in the scoring
                directly_associated_goterms = list(goterms_product_general)  # calling list constructor creates two separate entities, which prevents infinite looping !
                for directly_associated_goterm in directly_associated_goterms:  # WARNING: don't iterate over goterms_product_general, since this list is being updated in the for loop !!
                    parent_goterms = self.reverse_lookup.obo_parser.get_parent_terms(directly_associated_goterm)  # indirectly associated goterms
                    goterms_product_general += parent_goterms  # expand goterms_product_general by the parent goterms
                # delete duplicate parents by converting a list to set
                goterms_product_general = set(goterms_product_general)
                num_goterms_product_general = len(goterms_product_general)  # calculate new num_goterms_product_general

            for direction in ["+", "-"]:
                # num_goterms_product_process = sum(1 for goterm in process_goterms_list if (any(goterm_process['direction'] == direction for goterm_process in goterm.processes) and (any(product_id in goterm.products for product_id in product.id_synonyms) or product.genename in goterm.products)))
                # the above line is a single-line implementation of the below nested for loops
                num_goterms_product_process = 0  # all GO Terms which are associated with the current 'process' and the current 'direction' of regulation and are also associated with the current gene (product)
                goterms_product_process = [] # used for the results dict to be readable by the user
                goterms_product_process_ids = []  # used for calculation of child terms
                for goterm in process_goterms_list:  # iterate through each GO Term associated with the current pathophysiological process
                    for goterm_process in goterm.processes:  # goterm.processes holds which pathophysiological processes (eg. {'process': "cancer", 'direction': "-"}) the GO Term is associated with (this is determined by the user in the input.txt file)
                        if goterm_process["direction"] == direction:
                            if product.genename in goterm.products:  # attemp genename search first
                                num_goterms_product_process += 1
                                goterms_product_process.append(f"{goterm.id}: {goterm.name}")
                                goterms_product_process_ids.append(goterm.id)
                                break
                            for product_id in product.id_synonyms:  # if genename is not found, also look into product.id_synonyms
                                if product_id in goterm.products:
                                    num_goterms_product_process += 1
                                    goterms_product_process.append(f"{goterm.id}: {goterm.name}")
                                    goterms_product_process_ids.append(goterm.id)
                                    break

                # find all go term children of num_goterms_product_process here
                num_indirect_children = 0
                if self.reverse_lookup.model_settings.include_indirect_annotations is True:
                    all_indirect_children = []  # list of ids
                    for id in goterms_product_process_ids: # goterms_products_process_ids are all goterms from input.txt that are involved in the specific direction of regulation of a process (eg. cancer+)
                        child_goterms = self.reverse_lookup.obo_parser.get_child_terms(id)
                        all_indirect_children += child_goterms
                    all_indirect_children = set(all_indirect_children)  # to remove duplicates
                    # append to goterms product process so the user can see the results in data.json
                    for child_id in all_indirect_children:
                        goterms_product_process.append(f"indirect: {child_id}")
                    # update num_goterms_product_process
                    num_goterms_product_process += len(all_indirect_children)
                    num_indirect_children += len(all_indirect_children)

                # first, compute the sum of direct annotations
                num_goterms_all_process = sum(
                    1
                    for goterm in process_goterms_list
                    if any(
                        goterm_process["direction"] == direction
                        for goterm_process in goterm.processes
                    )
                )  # all of the GO Terms from input.txt file associated with the current process (and the process' regulation direction)
                num_goterms_all_process += num_indirect_children # if children are computed for goterms, then also increase num_goterms_all_process, to prevent negative values in the contingency table (specifically upper-right quadrant: num_goterms_all_process-num_goterms_product_process) 

                # time for Binomial test and "risk ratio"
                cont_table = [
                    [
                        num_goterms_product_process,
                        num_goterms_all_process - num_goterms_product_process
                    ],
                    [
                        num_goterms_product_general - num_goterms_product_process,
                        num_goterms_all_general - (num_goterms_all_process - num_goterms_product_process)
                    ],
                ]

                # check that any contingency table element is non-negative
                should_continue_current_loop = False
                for x in cont_table:
                    for y in x:
                        if y < 0:
                            stat_error = (
                                "Element of contingency table in class"
                                " fisher_exact_test is negative. All elements must be"
                                f" non-negative. Contingency table: {cont_table}. This"
                                " might be because a gene_name, which belongs to a"
                                " certain GO Term (obtained via web-download), isn't"
                                " found in the GO Annotations File."
                            )
                            results_dict[f"{process['process']}{direction}"] = {
                                # "n_prod_process" : num_goterms_product_process,
                                # "n_all_process" : num_goterms_all_process,
                                # "n_prod_general" : num_goterms_product_general,
                                # "n_all_general" : num_goterms_all_general,
                                "error": f"{stat_error}",
                                "num_terms_product_process": (
                                    num_goterms_product_process
                                ),
                                "num_terms_all_process": num_goterms_all_process,
                                "num_terms_product_general": (
                                    num_goterms_product_general
                                ),
                                "fold_enrichment": None,
                                "pvalue": None,
                                "odds_ratio": None,
                            }
                            should_continue_current_loop = True
                if (
                    should_continue_current_loop
                ):  # TODO: check the logic of should_continue_current_loop
                    continue

                fisher = fisher_exact(cont_table)
                fisher_pvalue = fisher.pvalue
                odds_ratio = fisher.statistic

                # calculate what amount of num_goterms_product_process would make this product statistically significant
                required_n_prod_process_for_stat_relevance = -1
                if D_DEBUG_CALCULATE_DESIRED_N_PROD_PROCESS and fisher_pvalue > 0.05:
                    previous_pvalue = fisher_pvalue
                    new_num_goterms_product_process = num_goterms_product_process
                    success = True  # if a pvalue less than 0.05 was found
                    while previous_pvalue > 0.05:
                        new_num_goterms_product_process += 1
                        new_cont_table = [
                            [
                                new_num_goterms_product_process,
                                num_goterms_all_process
                                - new_num_goterms_product_process,
                            ],
                            [
                                num_goterms_product_general
                                - new_num_goterms_product_process,
                                num_goterms_all_general
                                - (
                                    num_goterms_all_process
                                    - new_num_goterms_product_process
                                ),
                            ],
                        ]
                        all_non_negative = all(
                            value >= 0 for row in new_cont_table for value in row
                        )
                        if all_non_negative is False:
                            # logger.warning("Values in cont. table are negative.")
                            success = False
                            break
                        new_fisher = fisher_exact(new_cont_table)
                        new_pvalue = new_fisher.pvalue
                        if new_pvalue > previous_pvalue:
                            #logger.warning(
                            #    "Newly calculated pvalue is greater! prev_pvalue ="
                            #    f" {previous_pvalue}, new_pvalue = {new_pvalue}"
                            #)
                            success = False
                            break
                        else:
                            previous_pvalue = new_pvalue
                    if success:
                        required_n_prod_process_for_stat_relevance = (
                            new_num_goterms_product_process
                        )

                # TODO: NaN is necessary for further calculations! Do a json postproccess. START FROM HERE.
                # if math.isnan(odds_ratio) or math.isnan(fisher_pvalue):
                #    fisher_pvalue = None
                #    odds_ratio = None

                fold_enrichment_score = 0
                if (
                    num_goterms_all_process != 0
                    and num_goterms_product_general != 0
                    and num_goterms_all_general != 0
                ):
                    fold_enrichment_score = num_goterms_product_process / (
                        num_goterms_all_process
                        * (num_goterms_product_general / num_goterms_all_general)
                    )

                results_dict[f"{process['process']}{direction}"] = {
                    "n_prod_process": num_goterms_product_process,
                    "n_all_process": num_goterms_all_process,
                    "n_prod_general": num_goterms_product_general,
                    "n_all_general": num_goterms_all_general,
                    "num": num_goterms_product_process,
                    "required_n_prod_process_for_statistical_relevance": (
                        required_n_prod_process_for_stat_relevance
                    ),
                    "expected": num_goterms_all_process
                    * (
                        num_goterms_product_general / num_goterms_all_general
                        if num_goterms_all_general != 0
                        else 0
                    ),
                    "fold_enrichment": (
                        fold_enrichment_score
                    ),  # BUGFIX: ZeroDivisionError
                    "pvalue": fisher_pvalue,
                    "odds_ratio": odds_ratio,
                    "goterms_prod_process": goterms_product_process,
                }

        # all_target_pvalues = [results_dict[f"{process['process']}{process['direction']}"]['pvalue'] for process in self.reverse_lookup.target_processes]
        # combined_p = combine_pvalues(all_target_pvalues)
        # results_dict["comb_fisher_pvalue"] = combined_p.pvalue
        # combined_rr = statistics.mean([results_dict[f"{process['process']}{process['direction']}"]['odds_ratio'] for process in self.reverse_lookup.target_processes])
        # results_dict["comb_odds_ratio"] = combined_rr

        return results_dict


class inhibited_products_id(Metrics):
    """
    An implementation of the Metrics interface to return a list of all product ids inhibited by a specific miRNA, if the binding strength
    between a product id and a specific miRNA is greater than (ReverseLookup).miRNA_overlap_threshold.

    WARNING: field 'miRNA_overlap_threshold' must be defined in an instance of the ReverseLookup model passed to this constructor.

    Parameters:
      - (ReverseLookup) model: an instance of ReverseLookup

    Algorithm:
      Create an empty list.
      For the input miRNA, loop over all miRNA-mRNA binding strengths (stored in (miRNA).mRNA_overlaps)
        If binding strength > miRNA_overlap_threshold: append product id to list
      Return a list of product ids
    """

    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "inhibited_products_id"
        self.treshold = (
            self.reverse_lookup.miRNA_overlap_treshold
        )  # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold

    def metric(self, mirna: miRNA) -> List[str]:
        """
        An implementation of the scoring algorithm for a specific miRNA instance. It loops over all miRNA-mRNA binding strengths in (miRNA).mRNA_overlaps
        and returns a list of mRNA product ids, whose binding strengths to this miRNA are greater than miRNA_overlap_threshold.
        """
        inhibited_product_ids = []
        for product_id, overlap in mirna.mRNA_overlaps.items():
            if overlap >= self.treshold:
                inhibited_product_ids.append(product_id)
        return inhibited_product_ids


class basic_mirna_score(Metrics):
    """
    Score calculated from adv_score and overlap

    Scoring algorithm:
        Initialise score = 0.0
        For each product_id and it's float overlap (binding strength) value in (miRNA).mRNA_overlaps
        [TODO]: explain why the miRNA score is decreased (a = -1) if it binds to a product with a good threshold?
        if miRNA binds to a product well, then it's score should be increased!

    WARNING: [TODO] need to resolve scoring issue
    """

    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "basic_score"
        self.treshold = (
            self.reverse_lookup.miRNA_overlap_treshold
        )  # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold

    def metric(self, mirna: miRNA) -> float:
        """
        An implementation of the scoring algorithm for a specific miRNA instance. [TODO] explain more after the scoring issue is solved
        """
        score = 0.0
        for product_id, overlap in mirna.mRNA_overlaps.items():
            product = next(
                (x for x in self.reverse_lookup.products if x.uniprot_id == product_id),
                None,
            )  #  this line of code is looking through a sequence of products and finding the first product whose uniprot_id matches the value of product_id. If such a product is found, it is assigned to the variable product; otherwise, product is set to None
            if (
                product is not None
            ):  # each miRNA can have many products in it's 'mRNA_overlaps' field, this is a check that we are only analysing the products, which are also present in (ReverseLookup).products
                if overlap >= self.treshold:  # inhibited
                    a = (
                        -1
                    )  # deduct the score, since high score indicates the products is favourable for our target processes
                else:
                    a = 1
                score += a * product.scores["adv_score"]
        return score


class miRDB60predictor:
    # TODO: reimplement
    def __init__(self):
        # set the filepath to the miRDB prediction result file
        self._filepath = "app/goreverselookup/data_files/miRNAdbs/miRDB_v6.0_prediction_result.txt.gz"
        # check if the file exists and download it if necessary
        self._check_file()

        # read the file into memory and decode the bytes to utf-8
        # create a 2D dictionary between mRNAids, miRNA_ids (cols, rows) and their match_strengths
        self.mRNA_miRNA_match_strengths = {}
        self._readlines = []

        # TODO: this file opening mechanism is slow. it can surely be sped up.
        with gzip.open(self._filepath, "rb") as read_content:
            # self._readlines = [line.decode("utf-8") for line in read_content.readlines()]
            for line in read_content.readlines():
                line = line.decode("utf-8")
                self._readlines.append(line)

                miRNAid, mRNAid, match_strength = line.strip().split("\t")
                if mRNAid not in self.mRNA_miRNA_match_strengths:
                    self.mRNA_miRNA_match_strengths[
                        mRNAid
                    ] = (
                        {}
                    )  # if first-level dict doesn't exist, second-level (miRNAid) will throw an error
                self.mRNA_miRNA_match_strengths[mRNAid][miRNAid] = float(match_strength)

        # log the first 10 lines of the file
        # logger.info(self._readlines[:10])

        """
        for line in self._readlines:
            miRNAid, mRNAid, match_strength = line.strip().split("\t")
            if mRNAid not in self.mRNA_miRNA_match_strengths:
                self.mRNA_miRNA_match_strengths[mRNAid] = {} # if first-level dict doesn't exist, second-level (miRNAid) will throw an error
            self.mRNA_miRNA_match_strengths[mRNAid][miRNAid] = float(match_strength)
        """

    def _check_file(self):
        _max_retries = 5
        _retry_delay = 2
        # create the directory where the file will be saved if it doesn't exist
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            # download the file from the miRDB website with retries
            url = "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
            for i in range(_max_retries):
                try:
                    response = requests.get(url)
                    response.raise_for_status()
                except requests.exceptions.HTTPError as e:
                    logger.warning(f"Error downloading file from {url}: {e}")
                    if i < _max_retries - 1:
                        logger.warning(f"Retrying in {_retry_delay} seconds...")
                        time.sleep(_retry_delay)
                    else:
                        raise Exception(
                            f"Failed to download file from {url} after"
                            f" {self._max_retries} attempts"
                        ) from e

    def predict_from_product(
        self, product: Product, threshold: float = 0.0
    ) -> Dict[str, float]:
        """
        Finds all miRNAs and their match strengths (hsa-miR-xxx, 72.2) from miRDB_readlines for mRNA_refseq (e.g. NM_xxxxx).

        :param product: Product object with refseq_nt_id attribute
        :param threshold: Minimum match strength to include in result_list
        :return: Dictionary containing miRNAs as keys and match strengths as values
        """
        if not product.refseq_nt_id:
            return None

        result_dict = {}

        """ # update: changed to a 2D dictionary approach, which is faster than line by line
        # Iterate over each line in the list of read lines
        for line in self._readlines:
            # Check if product.refseq_nt_id is present in the line
            if product.refseq_nt_id.split(".")[0] in line:
                # Split the line by tabs to extract miRNA and match_strength
                miRNA, _, match_strength = line.strip().split("\t")
                # Convert match_strength to float
                match_strength = float(match_strength)
                # Add miRNA and match_strength to result_dict if match_strength >= threshold
                if match_strength >= threshold:
                    result_dict[miRNA] = match_strength
        """

        # 2D dictionary approach
        product_mRNA_id = product.refseq_nt_id.split(".")[0]
        if product_mRNA_id in self.mRNA_miRNA_match_strengths:
            corresponding_miRNA_matches = self.mRNA_miRNA_match_strengths[
                product_mRNA_id
            ]  # get dict of current mRNA strengths to different miRNAs, eg. {'hsa-miR-10393-3p': 60.3227953945, 'hsa-miR-122b-5p': 80.4854, 'hsa-miR-128-3p': 86.766, 'hsa-miR-1291': 64.6408466712, ...}
            for miRNA, match_strength in corresponding_miRNA_matches.items():
                if match_strength >= threshold:
                    result_dict[miRNA] = match_strength

        return result_dict

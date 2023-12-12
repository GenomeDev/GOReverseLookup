import os

from .Model import ReverseLookup
from .core.Metrics import (
    Metrics,
    adv_product_score,
    nterms,
    fisher_exact_test,
    binomial_test,
    basic_mirna_score,
)
from .core.Report import ReportGenerator
from .core.ModelSettings import ModelSettings

import logging
logger = logging.getLogger(__name__)
# from goreverselookup import logger


class WorkflowChecker:
    def __init__(self):
        self.steps = {
            "fetch_term_names": False,
            "fetch_term_products": False,
            "create_products": False,
            "fetch_uniprotid": False,
            "prune_products": False,
            "fetch_uniprot_infos": False,
            "score_products": False,
            "fetch_mRNA": False,
            "predict_miRNA": False,
            "score_miRNA": False,
        }

    def mark_step_done(self, step_name):
        if step_name in self.steps:
            self.steps[step_name] = True
        else:
            raise ValueError(f"Invalid step name: {step_name}")

    def is_step_done(self, step_name):
        if step_name in self.steps:
            return self.steps[step_name]
        else:
            raise ValueError(f"Invalid step name: {step_name}")


class Workflow:
    def __init__(
        self,
        input_file_fpath: str,
        save_folder_dir: str,
        model: ReverseLookup = None,
        name: str = "",
    ):
        """
        Creates a new Workflow instance.

        Parameters:
          - (str) input_file_fpath: The filepath to either the input.txt file or a data.json file for the model.
                                    If a 'model' is supplied, it will take precedence over creation of a model from input_file_fpath.
          - (str) save_folder_dir: The filepath to the directory, where the files will be saved
          - (ReverseLookup, optional) model: A ReverseLookup model. If no model is supplied, a new model will be created from input_file_fpath
          - (str) name: The name of this workflow
        """
        if model is None:
            if ".txt" in input_file_fpath:
                self.model = ReverseLookup.from_input_file(input_file_fpath)
            elif ".json" in input_file_fpath:
                self.model = ReverseLookup.load_model(input_file_fpath)
        else:
            self.model = model

        self.execution_sequence = []  # a list of functions to execute in sequence
        self.computed_scores = (
            {}
        )  # a dictionary between Metrics: (Metrics) aka metrics class - metrics instance of computed scores, computed scores are saved here from self.scores and cannot be deleted.
        self.scores = []  # a list of scoring algorithms, temporary, can be deleted
        self.goaf = self.model.goaf

        self.input_file_fpath = input_file_fpath
        self.save_folder_dir = save_folder_dir
        self.model_save_filepath = os.path.join(save_folder_dir, "data.json")
        self.model_statistically_relevant_products_filepath = os.path.join(save_folder_dir, "statistically_relevant_genes.json").replace("\\", "/")
        self.name = name

    def create_workflow():
        raise NotImplementedError(
            "You shouldn't call create_workflow on superclass Workflow. Create a"
            " subclass implementation of Workflow instead. using class"
            " WORKFLOW_NAME(Workflow)"
        )

    def add_function(self, function, *args, **kwargs):
        """
        Adds a function to self.execution_sequence.

        Usage:
        workflow = # instantiate a subclass Workflow instance, eg. WorkflowOne
        workflow.add_function(model.fetch_all_go_term_names_descriptions)
        workflow.add_function(model.fetch_all_go_term_products, web_download=True)
        workflow.run_workflow()
        """
        self.execution_sequence.append((function, args, kwargs))

    def perform_scoring(
        self,
        scoring_classes: Metrics | list[Metrics],
        delete_previous: bool = True,
        recalculate: bool = True,
    ):
        """
        This function instantiates each scoring class from 'scoring_classes' on the state of the current model instance (self.model).
        After all scoring classes are instantiated and stored into self.scores, model.score_products[self.scores] is called, to score
        all products in the current model instance with the supplied scores.

        Params:
          - (Metrics or list[Metrics]) scoring_classes: a single scoring class (an implementation of the Metrics class) or a list of
                                                        Metrics scoring algorithm implementations. Can be 'adv_product_score', 'nterms',
                                                        'fisher_exact_test', 'binomial_test', 'basic_mirna_score'
          - (bool) delete_previous: If True, will delete scoring results saved in self.scores. If False, will not delete scoring results
                                    saved in self.scores.
          - (bool) recalculate: If True, will recalculate scoring results, even if scoring results for product already exist. If False, won't
                                perform recalculations.

        Warning: When this function is re-called, self.scores will be deleted if delete_previous is True. All the computed scores are
        stored into the self.computed_scores dictionary, which maps relationships between a scoring algorithm's class to the latest
        instance of that class (aka the latest computation of the scores).

        Usage:
            # inside WorkflowOne -> create_workflow:
            self.perform_scoring([adv_score, nterms_score])
        """
        if delete_previous is True:
            self.scores = []  # clear previous scores

        if not isinstance(scoring_classes, list):
            scoring_classes = [scoring_classes]

        # instantiate all Metrics scoring algorithms on the current model state
        for scoring_class in scoring_classes:
            scoring_instance = None
            if isinstance(scoring_class, fisher_exact_test) or isinstance(
                scoring_class, binomial_test
            ):
                scoring_instance = scoring_class(self.model, self.goaf)
            else:  # if adv_product_score or nterms
                try:
                    scoring_instance = scoring_class(self.model)
                except TypeError as e:
                    # bugfix: binomial_test and fisher_exact are thrown into this clause somehow; maybe naming confusion error? TODO
                    if "missing" in str(e) and "positional argument" in str(e):
                        # missing positional argument -> supply goaf
                        scoring_instance = scoring_class(self.model, self.goaf)
                    else:
                        # raise other TypeErrors
                        raise e

            self.scores.append(scoring_instance)
            self.computed_scores[scoring_class] = scoring_instance

        # perform model scoring
        # note: if miRNA scoring, such as 'basic_miRNA_score' is used in self.scores, the model.score_products will redirect to model.score_miRNAs
        self.model.score_products(self.scores, recalculate=recalculate)

    def generate_report(
        self,
        product_scoring_algorithm: Metrics = None,
        miRNA_scoring_algorithm: Metrics = None,
        report_filepath: str = "",
        verbosity: int = 3,
    ):
        """
        Instantiates a Report object on the current state of self.model and generates the report.

        Parameters:
          - (Metrics, optional) product_scoring_algorithm: The scoring algorithm, the scores of which shall be used for the ordering of the products in the report.
                                                           If no product_scoring_algorithm is supplied, then this function will take the first element of self.computed_scores.
                                                           If self.comptued_scores is empty, a report will be generated without product scores.
          - (Metrics, optional) miRNA_scoring_algorithm: The scoring algorithm for miRNAs, which shall be used for the ordering of miRNAs in the report.
                                                         If no miRNA_scoring_algorithm is supplied, then this function will take the miRNA scoring algorithm of instance basic_mirna_score
                                                         from self.computed_scores. If it doesn't exist, then a report will be generated without miRNA scores.
          - (str) report_filepath: The filepath to the report .txt file, if left empty will default to self.save_folder_dir/report.txt
          - (int) verbosity: The verbosity of the report.
        """

        def pick_scoring_algorithm(
            scoring_algorithm_class: Metrics, miRNA: bool = False
        ):
            """
            Picks a scoring algorithm from self.computed_scores, which has the same instance as 'scoring_algorithm_class'.
            If 'scoring_algorithm_class' is None, then this function picks the first scoring algorithm it finds in self.computed_scores,
            irrespective of its type.

            Params:
              - (Metrics) scoring_algorithm_class: A class of subtype Metrics (eg. adv_product_score, nterms, fisher_exact_test, binomial_test, basic_mirna_score)
              - (bool) miRNA: if True, then it will attempt to return an instance of basic_miRNA_score
            """
            if scoring_algorithm_class is not None:
                for key, value in self.computed_scores.items():
                    if key == scoring_algorithm_class:
                        return value
                return None  # no key in self.comptued_scores matches scoring_algorithm_class = this score wasn't computed
            else:  # if scoring_algorithm_class is None
                for key, value in self.computed_scores.items():
                    if miRNA is True and isinstance(key, basic_mirna_score):
                        return value
                    elif miRNA is False:
                        return value
                return None

        final_report_filepath = ""
        if report_filepath == "":
            final_report_filepath = os.path.join(self.save_folder_dir, "report.txt")
        else:
            final_report_filepath = report_filepath

        # determine the product scoring algorithm and miRNA scoring algorithm based on the supplied classes
        final_product_scoring_algorithm = pick_scoring_algorithm(
            product_scoring_algorithm
        )
        final_miRNA_scoring_algorithm = pick_scoring_algorithm(
            miRNA_scoring_algorithm, miRNA=True
        )

        # generate the report
        report_generator = ReportGenerator(self.model, verbosity=verbosity)
        report_generator.general_report(
            final_report_filepath,
            final_product_scoring_algorithm,
            final_miRNA_scoring_algorithm,
        )

    def run_workflow(self):
        """
        Sequentially runs the functions specified in self.execution_sequence.
        """
        for function, args, kwargs in self.execution_sequence:
            function(*args, **kwargs)


class WorkflowOne(Workflow):
    def __init__(
        self,
        input_file_fpath: str,
        save_folder_dir: str,
        model: ReverseLookup = None,
        name: str = "",
        debug: bool = False,
    ):
        super().__init__(input_file_fpath, save_folder_dir, model, name)
        self.create_worklfow(debug=debug)

    def create_worklfow(self, debug: bool = False):
        # Pull GO data #
        # note: model = ReverseLookup.load_model(...) or ReverseLookup.load_from_file(...) is already handled in the Workflow constructor
        if debug:  # TODO: remove debug
            self.add_function(
                self.model._debug_shorten_GO_terms, 5
            )  # in debug mode, shorten the amount of GO Terms to 5

        self.add_function(self.model.fetch_all_go_term_names_descriptions)
        self.add_function(
            self.model.fetch_all_go_term_products,
            web_download=True,
            run_async=False,
            recalculate=False,
        )
        self.add_function(self.model.create_products_from_goterms)
        self.add_function(self.model.fetch_ortholog_products, refetch=False)
        self.add_function(self.model.save_model, self.model_save_filepath)
        self.add_function(self.model.fetch_product_infos, refetch=False)
        self.add_function(self.model.save_model, self.model_save_filepath)
        self.add_function(self.model.prune_products)
        self.add_function(self.model.save_model, self.model_save_filepath)

        # Score GO Term products #
        # self.add_function(self.perform_scoring, adv_product_score, nterms, fisher_exact_test, binomial_test)
        self.add_function(
            self.perform_scoring,
            scoring_classes=[adv_product_score, nterms],
            recalculate=False,
        )
        self.add_function(self.model.save_model, self.model_save_filepath)

        # Pull mRNA, perform mRNA-miRNA scoring #
        self.add_function(self.model.fetch_mRNA_sequences)
        self.add_function(self.model.save_model, self.model_save_filepath)
        self.add_function(self.model.predict_miRNAs)
        self.add_function(self.model.change_miRNA_overlap_treshold, 0.6, True)
        self.add_function(
            self.perform_scoring, scoring_classes=[basic_mirna_score], recalculate=False
        )
        self.add_function(self.model.save_model, self.model_save_filepath)

        # Generate report #
        self.add_function(
            self.generate_report,
            product_scoring_algorithm=adv_product_score,
            miRNA_scoring_algorithm=basic_mirna_score,
        )


class PrimaryWorkflow(Workflow):
    def __init__(
        self,
        input_file_fpath: str = "",
        save_folder_dir: str = "",
        model: ReverseLookup = None,
        name: str = "",
        debug: bool = False,
    ):
        # constructor chooses appropriate method to initialise the Model based on supplied parameters. A ReverseLookup 'model' instance takes precedence over input_file_fpath.
        super().__init__(input_file_fpath, save_folder_dir, model, name)
        self.create_workflow(debug=debug)

    def create_workflow(self, debug: bool = False, fetch_mirna=False):
        # Fetch all GO term names and descriptions
        self.add_function(
            self.model.fetch_all_go_term_names_descriptions, 
            run_async=True,
            req_delay=0.5,
            max_connections=20
        )
        # Fetch all GO term products
        self.add_function(
            self.model.fetch_all_go_term_products,
            web_download=True,
            run_async=True,
            recalculate=False,
            max_connections=40,
            request_params={"rows": 50000},
            delay=0.5,
        )
        # Create product instances from GO terms
        self.add_function(self.model.create_products_from_goterms)
        # Fetch human ortholog for products (either UniProtID, ENSG or genename)
        self.add_function(
            self.model.fetch_ortholog_products,
            refetch=False,
            run_async=True,
            max_connections=10,
            req_delay=0.1,
            semaphore_connections=5,
        )
        self.add_function(self.model.prune_products)
        self.add_function(self.model.save_model, self.model_save_filepath)
        # Fetch product information (from UniprotAPI or EnsemblAPI)
        self.add_function(
            self.model.fetch_product_infos,
            refetch=False,
            run_async=True,
            max_connections=15,
            semaphore_connections=10,
            req_delay=0.1,
        )
        self.add_function(self.model.prune_products)
        self.add_function(self.model.save_model, self.model_save_filepath)

        # Score products with the scores supplied in scoring_classes
        self.add_function(
            self.perform_scoring,
            scoring_classes=[
                adv_product_score,
                nterms,
                binomial_test,
                fisher_exact_test,
            ],
        )
        self.add_function(self.model.save_model, self.model_save_filepath)

        # Pull mRNA, perform mRNA-miRNA scoring
        if fetch_mirna is True:
            # Fetch mRNA sequences
            self.add_function(self.model.fetch_mRNA_sequences)
            self.add_function(self.model.save_model, self.model_save_filepath)
            # Predict miRNAs
            self.add_function(self.model.predict_miRNAs)
            # Score miRNAs
            self.add_function(self.model.change_miRNA_overlap_treshold, 0.6, True)
            self.add_function(
                self.perform_scoring,
                scoring_classes=[basic_mirna_score],
                recalculate=False,
            )
            self.add_function(self.model.save_model, self.model_save_filepath)

        # Perform statistical analysis of relevant products according to the chosen statistical test score
        # TODO: allow the user to choose which statistical test to use
        self.add_function(
            self.model.perform_statistical_analysis,
            test_name="fisher_test",
            filepath=self.model_statistically_relevant_products_filepath,
        )

        # TODO: generate report

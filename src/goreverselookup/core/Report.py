from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import os
from tabulate import tabulate
import traceback

if TYPE_CHECKING:
    from ..Model import ReverseLookup
    from .Metrics import Metrics
from ..util.FileUtil import FileUtil

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class ReportGenerator:
    def __init__(
        self,
        reverse_lookup: ReverseLookup,
        verbosity: int = 1,
        top_n: int = 5,
        width: int = 80,
    ):
        # initialize the report generator with a reverse lookup object and parameters for verbosity, top_n, and width
        self.reverse_lookup = (
            reverse_lookup  # the reverse lookup object to use for generating the report
        )
        self.width = width  # the width of the report output
        self.top_n = top_n  # the number of top results to display in the report
        self.verbosity = verbosity  # the level of detail to include in the report

    def _generate_header(self) -> str:
        """
        Generates the header part of the Report.py, currently hardcoded.
        """
        header = "Gene Ontology Reverse Lookup Tool".center(self.width) + "\n"
        header += (
            "Authors: Vladimir Smrkolj (SI), Aljosa Skorjanc (SI)".center(self.width)
            + "\n"
        )
        header += "March 2023".center(self.width) + "\n"
        return header

    def _generate_section(self, text: str) -> str:
        """
        Generates a specific section.
        Parameters:
          - (str) text: The section's heading text

        Returns:
          - (str) string: The formatted section heading text as displayed in the example.

        Example: The following is generated in the report:
        --------------------------------------------------------------------------------
                                    {text}
        --------------------------------------------------------------------------------
        """
        string = "-" * self.width + "\n"
        string += text.center(self.width) + "\n"
        string += "-" * self.width + "\n"
        return string

    def _generate_goterm_per_SOI_table(self) -> str:
        """
        This function creates a general overview table, crossing each SOI against positive (+), negative (-) or general (0) regulation,
        along with the count of GO Terms under each category. A second table is created, which again crosses each SOI against pos., neg. and
        general regulation, but included the GO IDs and GO labels of specific GO Terms under each category.

        Parameters: None

        Returns:
          - (str) string: The formatted tables using the data from the ReverseLookup model supplied in the constructor of this class

        Example:
        GO TERMS PER SOI
        +-----------+-----+-----+-----+---------+
        | SOI       |   + |   - |   0 |   Total |
        +===========+=====+=====+=====+=========+
        | diabetes  |   6 |   6 |  17 |      29 |
        +-----------+-----+-----+-----+---------+
        | angio     |   6 |   9 |  23 |      38 |
        +-----------+-----+-----+-----+---------+
        | Total     |  12 |  15 |  40 |      67 |
        +-----------+-----+-----+-----+---------+
        +-----------+--------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
        | SOI       | +                                                                                                                  | -                                                                                                                  | 0                                                                                                         |
        +===========+====================================================================================================================+====================================================================================================================+===========================================================================================================+
        | diabetes  | GO:1900077 - negative regulation of cellular response to insulin stimulus                                          | GO:1900078 - positive regulation of cellular response to insulin stimulus                                          | GO:0046950 - cellular ketone body metabolic process                                                       |
        |           | GO:0046627 - negative regulation of insulin receptor signaling pathway                                             | GO:0046628 - positive regulation of insulin receptor signaling pathway                                             | GO:0032869 - cellular response to insulin stimulus                                                        |
        |           | GO:0046676 - negative regulation of insulin secretion                                                              | GO:0032024 - positive regulation of insulin secretion                                                              | GO:0044381 - glucose import in response to insulin stimulus                                               |
        |           | GO:0043569 - negative regulation of insulin-like growth factor receptor signaling pathway                          | GO:0043568 - positive regulation of insulin-like growth factor receptor signaling pathway                          | GO:0006006 - glucose metabolic process                                                                    |
        |           | GO:0061179 - negative regulation of insulin secretion involved in cellular response to glucose stimulus            | GO:0038015 - positive regulation of insulin receptor signaling pathway by insulin receptor internalization         | GO:0043559 - insulin binding                                                                              |
        |           | GO:0038014 - negative regulation of insulin receptor signaling pathway by insulin receptor internalization         | GO:0035774 - positive regulation of insulin secretion involved in cellular response to glucose stimulus            | GO:1901143 - insulin catabolic process                                                                    |
        |           |                                                                                                                    |                                                                                                                    | GO:1901142 - insulin metabolic process                                                                    |
        |           |                                                                                                                    |                                                                                                                    | GO:0005899 - insulin receptor complex                                                                     |
        / ... /
        +-----------+--------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
        | angio     | GO:1903589 - positive regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis | GO:1903588 - negative regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis | GO:0001525 - angiogenesis                                                                                 |
        |           | GO:1903672 - positive regulation of sprouting angiogenesis                                                         | GO:1903671 - negative regulation of sprouting angiogenesis                                                         | GO:1903587 - regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis |
        |           | GO:0045766 - positive regulation of angiogenesis                                                                   | GO:0016525 - negative regulation of angiogenesis                                                                   | GO:1903670 - regulation of sprouting angiogenesis                                                         |
        |           | GO:0106090 - positive regulation of cell adhesion involved in sprouting angiogenesis                               | GO:0106089 - negative regulation of cell adhesion involved in sprouting angiogenesis                               | GO:0045765 - regulation of angiogenesis                                                                   |
        |           | GO:0090050 - positive regulation of cell migration involved in sprouting angiogenesis                              | GO:0090051 - negative regulation of cell migration involved in sprouting angiogenesis                              | GO:0060055 - angiogenesis involved in wound healing                                                       |
        |           | GO:0001935 - endothelial cell proliferation                                                                        | GO:0043532 - angiostatin binding                                                                                   | GO:0002040 - sprouting angiogenesis                                                                       |
        |           |                                                                                                                    | GO:0032311 - angiogenin-PRI complex                                                                                | GO:0002041 - intussusceptive angiogenesis                                                                 |
        / ... /
        +-----------+--------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
        """
        # Use a dictionary comprehension for the grouped_goterms initialization
        grouped_goterms = {
            (target["SOI"], direction): []
            for target in self.reverse_lookup.target_SOIs
            for direction in ("+", "-", "0")
        }

        # Use a for-loop to populate the grouped_goterms dictionary
        for goterm in self.reverse_lookup.goterms:
            for goterm_SOI in goterm.SOIs:
                key = (goterm_SOI["SOI"], goterm_SOI["direction"])
                grouped_goterms[key].append(goterm)

        string = "GO TERMS PER SOI" + "\n"

        if self.verbosity >= 1:
            table = [["SOI", "+", "-", "0", "Total"]]
            table.extend(
                [
                    [
                        target["SOI"],
                        len(grouped_goterms[(target["SOI"], "+")]),
                        len(grouped_goterms[(target["SOI"], "-")]),
                        len(grouped_goterms[(target["SOI"], "0")]),
                        sum(
                            [
                                len(grouped_goterms[(target["SOI"], "+")]),
                                len(grouped_goterms[(target["SOI"], "-")]),
                                len(grouped_goterms[(target["SOI"], "0")]),
                            ]
                        ),
                    ]
                    for target in self.reverse_lookup.target_SOIs
                ]
            )
            table.append(
                [
                    "Total",
                    sum([a[1] for a in table[1:]]),
                    sum([a[2] for a in table[1:]]),
                    sum([a[3] for a in table[1:]]),
                    sum([a[4] for a in table[1:]]),
                ]
            )
            string += (
                tabulate(table, headers="firstrow", tablefmt="grid").center(self.width)
                + "\n"
            )

        if self.verbosity == 2:
            table = [["SOI", "+", "-", "0"]]
            table.extend(
                [
                    [
                        target["SOI"],
                        "\n".join(
                            str(g.id) for g in grouped_goterms[(target["SOI"], "+")]
                        ),
                        "\n".join(
                            str(g.id) for g in grouped_goterms[(target["SOI"], "-")]
                        ),
                        "\n".join(
                            str(g.id) for g in grouped_goterms[(target["SOI"], "0")]
                        ),
                    ]
                    for target in self.reverse_lookup.target_SOIs
                ]
            )
            string += (
                tabulate(table, headers="firstrow", tablefmt="grid").center(self.width)
                + "\n"
            )

        if self.verbosity == 3:
            table = [["SOI", "+", "-", "0"]]
            table.extend(
                [
                    [
                        target["SOI"],
                        "\n".join(
                            str(f"{g.id} - {g.name}")
                            for g in grouped_goterms[(target["SOI"], "+")]
                        ),
                        "\n".join(
                            str(f"{g.id} - {g.name}")
                            for g in grouped_goterms[(target["SOI"], "-")]
                        ),
                        "\n".join(
                            str(f"{g.id} - {g.name}")
                            for g in grouped_goterms[(target["SOI"], "0")]
                        ),
                    ]
                    for target in self.reverse_lookup.target_SOIs
                ]
            )
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n"

        return string

    def _generate_goterms_statistics(self) -> str:
        """
        Generates general statistics for the GO Terms of the ReverseLookup model, supplied to the constructor of this class.

        Parameters:
          - None

        Returns:
          - (str) string: overall statistics for the model's GO Terms

        Example:
            GO TERMS STATISTICS
            Products per GO Term (min - avg - max): 0 - 63 - 501
        """

        string = "GO TERMS STATISTICS\n"

        # Calculate statistics and format the string
        if self.verbosity >= 1:
            products_per_goterm = [len(g.products) for g in self.reverse_lookup.goterms]
            min_g, max_g, avg_g = (
                min(products_per_goterm),
                max(products_per_goterm),
                sum(products_per_goterm) / len(products_per_goterm),
            )
            string += (
                f"Products per GO Term (min - avg - max): {min_g} - {avg_g:.0f} -"
                f" {max_g}\n"
            )

        return string

    def _generate_top_bottom_products_summary(self, score_key) -> str:
        """
        Generates (ReverseLookup).top_n (preset to 5) top-scored and bottom-scored products. The result is an overview table, where each of the
        top and bottom products are presented with their Gene Name, Score and Description. Then, (ReverseLookup).top_n * 2 tabels are also generated,
        one table per each presented top or bottom product, with all of it's GO Terms listed below, the GO Term Labels and GO Term Descriptions.

        The score_key corresponds to a scoring algorithm's name and is obtained from the product_score parameter supplied to the 'generate_report'
        function of this class, which should be one of the Metrics scoring implementations that was used to score the model.
        [TODO]: implement this more cleanly, this can be also pulled from the model's products after they are scored
        (eg. score_key = self.reverse_lookup.products[0].scores <- get the key of the scores dict !!)

        Returns:
          - (str) string: formatted result of top and bottom products

        Example:
        TOP and BOTTOM 5 PRODUCTS
        +-------------+---------+-------------------------------------------------------------+
        | Gene Name   | Score   | Description                                                 |
        +=============+=========+=============================================================+
        | LAMTOR4     | 18.20   | late endosomal/lysosomal adaptor, MAPK and MTOR activator 4 |
        +-------------+---------+-------------------------------------------------------------+
        | PROK1       | 18.20   | prokineticin 1                                              |
        +-------------+---------+-------------------------------------------------------------+
        | JCAD        | 15.46   | junctional cadherin 5 associated                            |
        +-------------+---------+-------------------------------------------------------------+
        | ACSBG1      | 15.40   | acyl-CoA synthetase bubblegum family member 1               |
        +-------------+---------+-------------------------------------------------------------+
        | TFDP2       | 15.40   | transcription factor Dp-2                                   |
        +-------------+---------+-------------------------------------------------------------+
        | ----        | ----    | ----                                                        |
        +-------------+---------+-------------------------------------------------------------+
        | GJA1        | -14.00  | gap junction protein alpha 1                                |
        +-------------+---------+-------------------------------------------------------------+
        | AGT         | -14.00  | angiotensinogen                                             |
        +-------------+---------+-------------------------------------------------------------+
        | LMNA        | -14.00  | lamin A/C                                                   |
        +-------------+---------+-------------------------------------------------------------+
        | C1QBP       | -15.40  | complement C1q binding protein                              |
        +-------------+---------+-------------------------------------------------------------+
        | ACAT2       | -15.40  | acetyl-CoA acetyltransferase 2                              |
        +-------------+---------+-------------------------------------------------------------+

                         LAMTOR4 - 18.20 - late endosomal/lysosomal adaptor, MAPK and MTOR activator 4
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        | GO term    | GO label                                                  | GO description                                                                                    |
        +============+===========================================================+===================================================================================================+
        | GO:0045766 | positive regulation of angiogenesis                       | Any process that activates or increases angiogenesis.                                             |
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        | GO:0046627 | negative regulation of insulin receptor signaling pathway | Any process that stops, prevents, or reduces the frequency, rate or extent of insulin receptor    |
        |            |                                                           | signaling.                                                                                        |
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        | GO:0045765 | regulation of angiogenesis                                | Any process that modulates the frequency, rate or extent of angiogenesis.                         |
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        | GO:0003158 | endothelium development                                   | The process whose specific outcome is the progression of an endothelium over time, from its       |
        |            |                                                           | formation to the mature structure. Endothelium refers to the layer of cells lining blood vessels, |
        |            |                                                           | lymphatics, the heart, and serous cavities, and is derived from bone marrow or mesoderm. Corneal  |
        |            |                                                           | endothelium is a special case, derived from neural crest cells.                                   |
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        | GO:0001938 | positive regulation of endothelial cell proliferation     | Any process that activates or increases the rate or extent of endothelial cell proliferation.     |
        +------------+-----------------------------------------------------------+---------------------------------------------------------------------------------------------------+
        / ... /
        """

        # Initialize the summary string with the header
        string = f"TOP and BOTTOM {self.top_n} PRODUCTS\n"

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(
            self.reverse_lookup.products,
            key=lambda x: x.scores[score_key],
            reverse=True,
        )
        top_products = sorted_products[: self.top_n]
        bottom_products = sorted_products[-self.top_n :]

        # If verbosity is at least 1, create a table of the top and bottom products with their scores and descriptions
        if self.verbosity >= 1:
            # Create the table as a list of lists and append each row
            table = [["Gene Name", "Score", "Description"]]
            for product in top_products:
                table.append(
                    [
                        product.genename,
                        f"{product.scores[score_key]:.2f}",
                        product.description,
                    ]
                )
            # Add a separator row and append each row for the bottom products
            table.append(["----", "----", "----"])
            for product in bottom_products:
                table.append(
                    [
                        product.genename,
                        f"{product.scores[score_key]:.2f}",
                        product.description,
                    ]
                )
            # Add the table to the summary string
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # If verbosity is at least 2, add details about the GO terms for each top and bottom product
        if self.verbosity >= 2:
            # Define a function to create the table of GO terms for a product
            def create_go_table(product):
                table_go = [["GO term", "GO label", "GO description"]]
                for goterm in self.reverse_lookup.get_all_goterms_for_product(product):
                    table_go.append([goterm.id, goterm.name, goterm.description])
                return table_go

            # Add the details for the top products
            for product in top_products:
                string += (
                    " " * 10 + f"{product.genename} - {product.scores[score_key]:.2f} -"
                    f" {product.description}".center(100) + "\n"
                )
                string += (
                    tabulate(
                        create_go_table(product),
                        headers="firstrow",
                        tablefmt="grid",
                        maxcolwidths=100,
                    )
                    + "\n\n"
                )

            # Add a separator and add the details for the bottom products
            string += ("-" * 30) + "\n\n"
            for product in bottom_products:
                string += (
                    " " * 10 + f"{product.genename} - {product.scores[score_key]:.2f} -"
                    f" {product.description}".center(100) + "\n"
                )
                string += (
                    tabulate(
                        create_go_table(product),
                        headers="firstrow",
                        tablefmt="grid",
                        maxcolwidths=100,
                    )
                    + "\n\n"
                )

        return string

    def _generate_top_miRNAs_summary(self, products_score_key, mirna_score_key) -> str:
        """
        TODO: comment
        """
        # Create the header string.
        string = f"TOP {self.top_n} miRNAs" + "\n"
        string += (
            f"+ annotates Product which is in top {self.top_n}, and - annotates Product"
            f" which is in bottom {self.top_n}\n\n"
        )

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(
            self.reverse_lookup.products,
            key=lambda x: x.scores[products_score_key],
            reverse=True,
        )
        top_products = sorted_products[: self.top_n]
        bottom_products = sorted_products[-self.top_n :]

        # Get the top miRNAs.
        top_miRNAs = sorted(
            self.reverse_lookup.miRNAs,
            key=lambda x: x.scores[mirna_score_key],
            reverse=True,
        )[: self.top_n]

        # If verbosity is set to 1, create a table with the top miRNAs and their scores.
        if self.verbosity == 1:
            table = [["miRNA", "mirna_score_key"]]
            for _miRNA in top_miRNAs:
                table.append([_miRNA.id, _miRNA.scores[mirna_score_key]])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # If verbosity is set to 2 or higher, create a table with the top miRNAs, their scores, and the products they inhibit.
        if self.verbosity >= 2:
            table = [["miRNA", "mirna_score_key", "suppressed products"]]

            for _miRNA in top_miRNAs:
                inhibited_product_id = []
                for product_id, overlap in _miRNA.mRNA_overlaps.items():
                    if overlap >= self.reverse_lookup.miRNA_overlap_treshold:
                        inhibited_product_id.append(product_id)

                temp_list = []
                temp_t_list = []
                temp_b_list = []
                for product_id in inhibited_product_id:
                    if (
                        product_id is None or product_id == "null"
                    ):  # TODO: CONTINUE FROM HERE! CHECK THAT THIS DOESN'T BREAK THE REPORT / CAUSE THE REPORT GENERATOR TO REPORT 0 SUPPRESSED ELEMENTS.
                        continue
                    if any(product_id in sub.uniprot_id for sub in top_products):
                        temp_t_list.append(f"+{product_id}")
                    if any(product_id in sub.uniprot_id for sub in bottom_products):
                        temp_b_list.append(f"-{product_id}")
                    else:
                        temp_list.append(f"{product_id}")
                temp_list = (
                    temp_t_list + temp_list
                )  # To nsure the top/bottom products are displyed on top.
                temp_list = temp_b_list + temp_list
                table.append(
                    [
                        _miRNA.id,
                        _miRNA.scores[mirna_score_key],
                        "\n".join(item for item in temp_list),
                    ]
                )
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # Return the summary string.
        return string

    def general_report(
        self,
        filepath: str,
        product_score: Optional[Metrics] = None,
        miRNA_score: Optional[Metrics] = None,
    ):
        """
        Generates the general report and writes it to a file.

        Args:
          - (str) filepath: The path to the output file.
          - (list(Metrics)) product_score: one or more Metrics implementations, which were used for the scoring of products
          - (list(Metrics)) miRNA_score: one ore more Metrics implementations, which were used for the scoring of miRNAs

        Example usage:
            # Pull GO data #
            model = ReverseLookup.load_model("diabetes_angio_2/data.json")
            model.fetch_all_go_term_names_descriptions()
            model.fetch_all_go_term_products(web_download=True)
            model.create_products_from_goterms()
            model.fetch_ortholog_products(refetch=False)
            model.fetch_product_infos(refetch=False)
            model.prune_products()
            model.save_model("diabetes_angio_2/data.json")

            # Score GO Term products #
            adv_score = adv_product_score(model)
            nterms_score = nterms(model)
            model.score_products([adv_score, nterms_score])
            model.save_model("diabetes_angio_2/data.json")

            # Score mRNA-miRNA #
            model.fetch_mRNA_sequences()
            model.predict_miRNAs()
            model.change_miRNA_overlap_treshold(0.6, True)
            basic_score = basic_mirna_score(model)
            model.score_miRNAs(basic_score)
            model.save_model("diabetes_angio_1/data.json")

            # Generate report #
            report = ReportGenerator(model, verbosity=3)
            report.general_report("diabetes_angio_1/general.txt", product_score=adv_score)
        """
        filepath = filepath.replace(
            "/", os.sep
        )  # # replace any '/' to avoid having both \\ and / in a single filepath

        # Generate header of the report
        report = self._generate_header() + "\n\n"

        # Generate section on GOTerms
        if len(self.reverse_lookup.goterms) > 0:
            report += self._generate_section("GO TERMS")
            report += self._generate_goterms_statistics() + "\n"
            report += self._generate_goterm_per_SOI_table() + "\n"

        # Generate section on Products
        if len(self.reverse_lookup.products) > 0:
            report += self._generate_section("PRODUCTS")
            # TODO: bottom line results in error, if no product_score is supplied !!!
            report += (
                self._generate_top_bottom_products_summary(product_score.name) + "\n"
            )

        # Generate section on miRNAs
        if len(self.reverse_lookup.miRNAs) > 0:
            report += self._generate_section("miRNAs")
            report += (
                self._generate_top_miRNAs_summary(product_score.name, miRNA_score.name)
                + "\n"
            )

        if not os.path.isabs(filepath):
            current_dir = os.path.dirname(os.path.abspath(traceback.extract_stack()[0].filename))
            mac_filepath = os.path.join(current_dir, filepath)  # mac_filepath, since this approach works on a mac computer

        # Create directory for the report file, if it does not exist
        try:
            os.makedirs(
                os.path.dirname(mac_filepath), exist_ok=True
            )  # this approach works on a mac computer

            # Write the report to the output file
            with open(mac_filepath, "w") as f:
                f.write(report)
        except OSError:
            # TODO
            # first pass is allowed, on Windows 10 this tries to create a file at
            # 'C:\\Program Files\\Python310\\lib\\diabetes_angio_1/general.txt'
            # which raises a permission error.
            pass

        # fallback if the above fails
        try:
            # current_dir = os.getcwd()
            # win_filepath = os.path.join(current_dir, filepath) # reassign filepath to absolute path
            win_filepath = FileUtil.find_win_abs_filepath(filepath)
            os.makedirs(os.path.dirname(win_filepath), exist_ok=True)

            with open(win_filepath, "w") as f:
                f.write(report)
        except OSError:
            logger.info(f"ERROR! Cannot make directory at {os.path.dirname(filepath)}")

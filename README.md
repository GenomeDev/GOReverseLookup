# GOReverseLookup
**GOReverseLookup** GOReverseLookup is a Python package designed for Gene Ontology Reverse Lookup. It serves the purpose of identifying statistically significant genes within a set of selected Gene Ontology Terms.

While Gene Ontology offers valuable insights through gene annotations associated with individual terms, the biological reality often involves the complex interaction of multiple terms that either promote or inhibit a specific pathophysiological process. Unfortunately, Gene Ontology does not provide a built-in mechanism for computing statistically significant genes that are shared across multiple terms.

GOReverseLookup steps in to bridge this gap. It empowers researchers to uncover genes of statistical significance that participate in various interconnected Gene Ontology Terms, shedding light on intricate biological processes.

# Getting started
This section instructs you how to install the GOReverseLookup package and it's prerequisites.
## Prerequisites
* Python >= 3.10.0
* Integrated Development Application (IDE) - such as Visual Studio Code (VSCode)
## Installation
To use the package, install it with pip:
```
pip install goreverselookup
```

# Usage
## Using an IDE (e.g., Visual Studio Code)
### Folder structure setup
Firstly, create a folder that will be the root of your Python project. We will refer to the root folder as `root/`. Ideally, you should create the following folder structure:
```
root/
    - input_files
        - input.txt
    - results
        - data.json
        - statistically_significant_genes.json
    - main.py
```
Explanation of the folder structure:
    - `input_files` is where input.txt files are stored. These files will server as the entry point for the program and will have to be constructed manually.
    - `results` is where the results of the analysis will be stored. After the program runs in entirety, two files will be computed: 
        - `data.json` 
        - `statistically_significant_genes.json`
        The contents of the above files are explained in subsequent sections.
    - `main.py` is the main file, where the Python code will be placed to carry out the analysis.

### Code
There are two main algorithms that can be used to achieve the same result. The mentioned scripts should be put in the `main.py` file, which you should execute using a IDE of your choice, preferably Visual Studio Code.

a) ** The Workflows Algorithm**
Workflows provide a simple and easy-to-use solution to jumpstart your research. All that is needed is that you provide an input file and a save folder to the PrimaryWorkflow class and call workflow.run_workflow():
```
# import necessary classes
from goreverselookup import PrimaryWorkflow
from goreverselookup import Cacher
from goreverselookup import LogConfigLoader

# setup logger
import logging
LogConfigLoader.setup_logging_config()
logger = logging.getLogger(__name__)

# setup cacher
Cacher.init(cache_dir="app/goreverselookup/cache")

# run the research algorithm
workflow = PrimaryWorkflow(input_file_fpath="input_files/input.txt", save_folder_dir="results")
workflow.run_workflow()
```
PrimaryWorkflow expects two parameters: the first is `input_file_fpath`, which should be set to the path of your input.txt file. The second parameter is `save_folder_dir`, which is a path to the directory where the result files will be saved.

If the workflow executes successfully, two files should be saved into `save_folder_dir`:
    - `data.json` contains the representation of the entire research workflow (model) and can be later used to load the model instead of recomputing it from input.txt again
    - `statistically_significant_genes.json` contains the genes which were found to statistically significantly contribute to the processes the researcher is interested in

**Logger** is set up in order to log the current algorithm steps to the console of your IDE (e.g., VSCode). By setting up the logger, you can monitor which commands the program is currently executing.

**Cacher** is set up in order to save web responses and and function return values. It's function is to speed up any subsequent program runs on data that has been already requested from the servers or data already computed in a previous time point. This is useful, because if you make slight modifications to input.txt and re-run the program, Cacher will restore the already computed values from the program's cache and only send the web requests only for the new additions to the input.txt file. 

b) ** The ReverseLookup Model Algorithm **
The ReverseLookup (research) model is the core of the analysis. Workflows are actually just a wrapper around the ReverseLookup model, hiding all the complexities from the researcher. To carry out the same full analysis (as is done in the workflows algorithm), you need to construct and run the following `main.py` file:
```
from goreverselookup import ReverseLookup
from goreverselookup import Cacher
from goreverselookup import LogConfigLoader
from goreverselookup import nterms, adv_product_score, binomial_test, fisher_exact_test

# setup logger
import logging
LogConfigLoader.setup_logging_config()
logger = logging.getLogger(__name__)

# setup cacher
Cacher.init(cache_dir="app/goreverselookup/cache")

# load the model from input file and query relevant data from the web
model = ReverseLookup.from_input_file("input_files/input.txt")
model.fetch_all_go_term_names_descriptions(run_async=True, req_delay=0.1)
model.fetch_all_go_term_products(web_download=True, run_async=True)
model.create_products_from_goterms()
model.fetch_ortholog_products(run_async=True, max_connections=15, semaphore_connections=5, req_delay=0.1)
model.prune_products()
model.save_model("results/data.json")

# test model load from existing json, perform model scoring
model = ReverseLookup.load_model("results/data.json")
nterms_score = nterms(model)
adv_prod_score = adv_product_score(model)
binom_score = binomial_test(model)
fisher_score = fisher_exact_test(model)
model.score_products(score_classes=[nterms_score, adv_prod_score, binom_score, fisher_score])
model.model_settings.pvalue = 0.05 # set pvalue to be used in statistical analysis
model.perform_statistical_analysis(test_name="fisher_test", filepath="results/statistically_relevant_genes.json")
model.save_model("results/data.json")
```
Each function call on the `ReverseLookup` instance called `model` has a descriptive name, highlighting it's task. The functions are heavily commented in code, so we encourage you to explore these comments and the code of our GitHub repository so as to gain a deeper understanding of the inner complexities of our tool.

## Using GOReverseLookup.exe
[todo]


# Roadmap
[todo]: link github projects here

# Contributing
Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

# License
Distributed using Apache 2.0 License. See `LICENSE.txt` for more information.

# Contact
Aljoša Škorjanc - skorjanc.aljosa@gmail.com
Vladimir Smrkolj - smrkolj.vladimir@gmail.com

# GOReverseLookup

[![PyPI package](https://img.shields.io/badge/pip%20install-goreverselookup-brightgreen)](https://pypi.org/project/goreverselookup/) [![version number](https://img.shields.io/github/v/release/MediWizards/GOReverseLookup)](https://github.com/MediWizards/GOReverseLookup/releases) [![Actions Status](https://img.shields.io/github/actions/workflow/status/MediWizards/GOReverseLookup/test_on_push.yml)](https://github.com/MediWizards/GOReverseLookup/actions/workflows/test_on_push.yml) [![License](https://img.shields.io/github/license/MediWizards/GOReverseLookup)](https://github.com/MediWizards/GOReverseLookup/blob/main/LICENSE)

**GOReverseLookup** is a Python package designed for Gene Ontology Reverse Lookup. It serves the purpose of identifying statistically significant genes within a set or a cross-section of selected Gene Ontology Terms. Researchers need only define their own states of interest (SOIs), and select GO terms must be attributed as either positive or negative regulators of the chosen SOIs. For more information regarding the creation of the input file for the program, refer to the _Input file_ section. Once the input file is created, the GOReverseLookup program can be started. Once the algorithm is completed, the program saves statistically significant genes in a standalone file.

For example, if researchers were interested in the _angiogenesis_ SOI, then an attributed group of GO terms as positive regulators of _angiogenesis_ might have been defined using the following GO terms:
- GO:1903672 positive regulation of sprouting angiogenesis
- GO:0001570 vasculogenesis
- GO:0035476 angioblast cell migration

And negative regulators of the _angiogenesis_ SOI might have been defined as the following group:
- GO:1903671 negative regulation of sprouting angiogenesis
- GO:1905554 negative regulation of vessel branching
- GO:0043537 negative regulation of blood vessel endothelial cell migration

If a researcher defines the _target process_ as positive regulation of a desired SOI (in our case _angiogenesis_), then GOReverseLookup finds all genes statistically relevant for the group of GO terms defined as positive regulators of _angiogenesis_ (p < 0.05) while excluding any genes determined to be statistically significant (_p_ < 0.05) in the opposing process (in our case, negative regulation of _angiogenesis_). _P_-value threshold can also be manually set by the user.


## Getting Started
This section instructs you how to install the GOReverseLookup package and its prerequisites.

### Folder setup
You MUST create a local folder anywhere on your disk, which will be used as the GOReverseLookup's working environment, as well as unified storage for all of your research projects. We advise you to create a folder structure with a folder named `goreverselookup` as the parent folder (this folder will be used as a local installation location for the GOReverseLookup program), and a subfolder named `research_models`, where you will store the input files for GOReverseLookup and their results. Therefore, the folder structure should be the following:
```
.../goreverselookup/
    - research_models/
```

### Installation
#### Python installation
For your computer to understand the GOReverseLookup program, it requires the Python programming language, which MUST be installed. Our program is currently tested on Python versions 3.10.x through 3.11.x, but not yet on 3.12.x. Thus, we advise you to use the Python version 3.11.5, which is available for download from [this website](https://www.python.org/downloads/release/python-3115/). Following this link, navigate to the _Files_ section:
- if you are using Windows: download _Windows installer (64-bit)_
- if you are using macOS: download _macOS 64-bit universal2 installer_

![github python Files section](https://i.ibb.co/kXLg2QD/goreverselookup-pyth.png)

Open the File Explorer program, then open the Downloads folder and run the installer by double clicking it.

![downloads folder](https://i.ibb.co/8xFzpjY/github-downloads-folder.png)
![python installer](https://i.ibb.co/JcnB96N/github-pyinstall.png)

The default Python installer window pops up:

<img src="https://i.ibb.co/YR8qZMc/github-pyinstaller-greet.png" width="450">

**Make sure** to also select **Add python.exe to PATH**. This will make Python available across all-file locations, which is of extreme importance for running Python commands from the console (Command prompt in Windows). Then, click on **Install Now**. A further observation of the installer's window also reveals that this installer is bundled with PIP (Python's package manager), thus manual installation of PIP won't be necessary. This is important, since PIP will be used to download GOReverseLookup.

<img src="https://i.ibb.co/5R9TrxP/github-pyinstaller-addpath.png" width="450">

Wait for the installation of Python to finish. Once it is finished, close the installer window.

If you wish to download a specific Python version, browse through the [Python's downloads page](https://www.python.org/downloads/) - for beginners, we advise you to find a release with an available installer. 

Then, open the command prompt using the Windows search bar:

<img src="https://i.ibb.co/PjkK65L/github-cmd.png" width="450">

Inside the command prompt, execute the command `python --version`. If Python installation has been completed successfully, a version of the Python programming language will be displayed:

![github cmd python version](https://i.ibb.co/RY55LKv/github-cmd-pyvers.png)

Also verify that PIP (Python's package manager) is installed. In our instance, it has been mentioned in the Python installer's window that PIP will also be installed along with Python. To verify the installation of PIP, run the `pip --version` command:

![cmd pip version](https://i.ibb.co/NFNgL40/github-cmd-pipvers.png)

#### Creating your GOReverseLookup workspace
To create a standalone GOReverseLookup workspace that will be central both to GOReverseLookup's installation files and the research files, create the folder setup as instructed in _Folder setup_. Create a Python's virtual environment in the `goreverselookup` folder using the command `python -m venv "PATH_TO_GOREVERSELOOKUP"`. For example, on my computer, the `goreverselookup` folder exists at `F:\Development\python_environments\goreverselookup`, thus the command to create the virtual environment is: `python -m venv "F:\Development\python_environments\goreverselookup"`:

![pyvenv](https://i.ibb.co/3pNDHF0/github-pyvenv.png)

To find the path to your goreverselookup folder, open the goreverselookup folder in the File Explorer and click on the Address Bar, then copy the filepath.

![goreverselookup file explorer path 1](https://i.ibb.co/f2qp6Lr/github-fe-goreverselookup.png)

![goreverselookup file explorer path 2](https://i.ibb.co/3dwZkyd/github-fe-adrbar.png)

After running the virtual environment creation command, you should notice the goreverselookup folder be populated with new folders: `Include`, `Lib` and `Scripts`, and a file named `pyvenv.cfg`. These belong to the newly created Python's virtual environment, so do not change their contents in any way. As stated in the _Folder setup_ section, the goreverselookup folder also contains a `research_models` folder.

![goreverselookup folder after pyvenv](https://i.ibb.co/D7MfZVv/github-grvfolder-after-pyvenv.png)

To activate the newly created virtual environment, there exists an activation script named `activate.bat` in the newly created `Scripts` folder. You will need to activate this virtual environment in command prompt every time you begin working with GOReverseLookup, thus we advise you to save the activation command in a text file somewhere easily accessible, such as your desktop. To activate the virtual environment, just supply the path to the activation script to the command prompt - in our case, the path to the activation script is `F:\Development\python_environments\goreverselookup\Scripts\activate`. After running this in command prompt, the virtual environment will be activated:

![goreverselookup venv activation](https://i.ibb.co/Gx7g8kF/github-venvactivation.png)

#### Installing GOReverseLookup
As per instructions in _Creating your GOReverseLookup workspace_, activate the newly created virtual environment, so the current command prompt pointer points to the virtual environment. E.g.:

<img src="https://i.ibb.co/Gx7g8kF/github-venvactivation.png" width="350">

Now, run the command `pip install goreverselookup` and wait for the installation to complete:

![goreverselookup pip install](https://i.ibb.co/T2drrGF/github-goreverselookup-install-pip.png)

To confirm the installation, run the command `pip list` and find the `goreverselookup` package, along with it's version:

![goreverselookup pip list](https://i.ibb.co/R99Rp9H/github-pip-list.png)


## Usage
### Creating the input file
The entry to the program is an input file, which is ideally placed in the `.../goreverselookup/research_models/` folder, as explained in _Folder setup_. It contains all the relevant data for the program to complete the analysis of statistically important genes that positively or negatively contribute to one or more states of interest.

An example input.txt file to discover the genes that positively contribute to both the development of chronic inflammation and cancer is:
```
# Comments are preceded by a single '#'. Comment lines will not be parsed in code.
# Section titles are preceded by three '###'
# The values at each line are usually delineated using the TAB character. E.g. pvalue    0.05 (pvalue and it's value 0.05 are separated by a TAB).
#
###evidence_code_groups
experimental	EXP_ECO:0000269,IDA_ECO:0000314,IPI_ECO:0000353,IMP_ECO:0000315,IGI_ECO:0000316,IEP_ECO:0000270,HTP_ECO:0006056,HDA_ECO:0007005,HMP_ECO:0007001,HGI_ECO:0007003,HEP_ECO:0007007
phylogenetic	IBA_ECO:0000318,IBD_ECO:0000319,IKR_ECO:0000320,IRD_ECO:0000321
computational_analysis	ISS_ECO:0000250,ISO_ECO:0000266,ISA_ECO:0000247,ISM_ECO:0000255,IGC_ECO:0000317,RCA_ECO:0000245
author_statement	TAS_ECO:0000304,NAS_ECO:0000303
curator_statement	IC_ECO:0000305,ND_ECO:0000307
electronic	IEA_ECO:0000501
###settings
pvalue	0.05
target_organism	homo_sapiens|UniProtKB|NCBITaxon:9606 # format: organism_label|organism_database|ncbi_taxon
ortholog_organisms	danio_rerio|ZFIN|NCBITaxon:7955,rattus_norvegicus|RGD|NCBITaxon:10116,mus_musculus|MGI|NCBITaxon:10090,xenopus_tropicalis|Xenbase|NCBITaxon:8364
evidence_codes	experimental(~),phylogenetic(~),computational_analysis(~),author_statement(TAS),!curator_statement(ND),!electronic(~)
#evidence_codes	experimental(~),phylogenetic(~),computational_analysis(~),author_statement(TAS),!curator_statement(ND),electronic(~)
gorth_ortholog_fetch_for_indefinitive_orthologs	True
gorth_ortholog_refetch	False
fisher_test_use_online_query	False
include_indirect_annotations	False
uniprotkb_genename_online_query	False
goterm_gene_query_timeout	240
goterm_gene_query_max_retries	3
exclude_opposite_regulation_direction_check	False
###filepaths
go_obo	data_files/go.obo	https://purl.obolibrary.org/obo/go.obo	all
goa_human	data_files/goa_human.gaf	http://geneontology.org/gene-associations/goa_human.gaf.gz	homo_sapiens
#goa_zfin TODO
#goa_rgd TODO
#goa_mgi TODO
#goa_xenbase TODO
ortho_mapping_zfin_human	data_files/zfin_human_ortholog_mapping.txt	https://zfin.org/downloads/human_orthos.txt	danio_rerio
ortho_mapping_mgi_human	data_files/mgi_human_ortholog_mapping.txt	https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt	mus_musculus
ortho_mapping_rgd_human	data_files/rgd_human_ortholog_mapping.txt	https://download.rgd.mcw.edu/data_release/HUMAN/ORTHOLOGS_HUMAN.txt	rattus_norvegicus
ortho_mapping_xenbase_human	data_files/xenbase_human_ortholog_mapping.txt	https://download.xenbase.org/xenbase/GenePageReports/XenbaseGeneHumanOrthologMapping.txt	xenopus
###processes [proces name] [to be expressed + or suppressed -]
chronic_inflammation	+
cancer	+
###categories [category] [True / False]
biological_process	True
molecular_activity	True
cellular_component	False
###GO_terms [GO id] [process] [upregulated + or downregulated - or general 0] [weight 0-1] [GO term name - optional] [GO term description - optional]
GO:0006954	chronic_inflammation	+	1	inflammatory response
GO:1900408	chronic_inflammation	-	1	negative regulation of cellular response to oxidative stress
GO:1900409	chronic_inflammation	+	1	positive regulation of cellular response to oxidative stress
GO:2000524	chronic_inflammation	-	1	negative regulation of T cell costimulation
GO:2000525	chronic_inflammation	+	1	positive regulation of T cell costimulation
GO:0002578	chronic_inflammation	-	1	negative regulation of antigen processing and presentation
GO:0002579	chronic_inflammation	+	1	positive regulation of antigen processing and presentation
GO:1900017	chronic_inflammation	+	1	positive regulation of cytokine production involved in inflammatory response
GO:1900016	chronic_inflammation	-	1	negative regulation of cytokine production involved in inflammatory response
GO:0001819	chronic_inflammation	+	1	positive regulation of cytokine production
GO:0001818	chronic_inflammation	-	1	negative regulation of cytokine production
GO:0050777	chronic_inflammation	-	1	negative regulation of immune response
GO:0050778	chronic_inflammation	+	1	positive regulation of immune response
GO:0002623	chronic_inflammation	-	1	negative regulation of B cell antigen processing and presentation
GO:0002624	chronic_inflammation	+	1	positive regulation of B cell antigen processing and presentation
GO:0002626	chronic_inflammation	-	1	negative regulation of T cell antigen processing and presentation
GO:0002627	chronic_inflammation	+	1	positive regulation of T cell antigen processing and presentation

GO:0007162	cancer	+	1	negative regulation of cell adhesion
GO:0045785	cancer	-	1	positive regulation of cell adhesion
GO:0010648	cancer	+	1	negative regulation of cell communication
GO:0010647	cancer	-	1	positive regulation of cell communication
GO:0045786	cancer	-	1	negative regulation of cell cycle
GO:0045787	cancer	+	1	positive regulation of cell cycle
GO:0051782	cancer	-	1	negative regulation of cell division
GO:0051781	cancer	+	1	positive regulation of cell division
GO:0030308	cancer	-	1	negative regulation of cell growth
GO:0030307	cancer	+	1	positive regulation of cell growth
#GO:0043065	cancer	-	1	positive regulation of apoptotic process
#GO:0043066	cancer	+	1	negative regulation of apoptotic process
GO:0008285	cancer	-	1	negative regulation of cell population proliferation
GO:0008284	cancer	+	1	positive regulation of cell population proliferation
```





### Dependencies

* Describe any prerequisites, libraries, OS version, etc., needed before installing program.
* ex. Windows 10

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)

Known limitations:
When using asynchronous querying for GO term products, if one of the requests inside a batch of requests exceeds the 'goterm_gene_query' timeout value (one of the settings), the entire batch of product queries will fail. This usually happens when the user attempts to collect products of GO terms with millions of more annotated genes. For us, an experimental 'goterm_gene_query' timeout value that successfully queris GO terms with ~1 million annotated genes is 240 seconds.

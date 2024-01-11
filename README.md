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
###states_of_interest [SOI name] [to be expressed + or suppressed -]
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
The main role of the researcher is to establish one or more custom states of interest (SOIs) and then attribute specific GO terms to the SOIs. Thus, SOIs and GO term attributions will be covered first.

#### Creating SOIs (states_of_interest section)
States of interest are created in the `states_of_interest` section. A SOI represents a name of a specific state of interest. Besides the name, either `+` or `-` is added in the line besides the SOI name in order to specify whether the researcher is interested to find genes responsible in the positive contribution (stimulation) of the SOI or negative contribution (inhibition) of the SOI.

Example: A researcher observes increased capillary growth in a histological sample. An example SOI could be `angiogenesis    +` to discover genes involved in the stimulation of angiogenesis.

#### Attributing GO terms to SOIs (GO_terms section)
After SOIs have been created, they need attributed GO terms to specifically define them. GO terms are attributed to SOIs in the `GO_terms` section, by first naming a GO term id, followed by the SOI, the impact of the GO term on the SOI (+ or -), a weight (this is historical and is kept at 1) and a description of the GO term.

Example: A researcher defined an `angiogenesis` SOI. Now, the researcher can assign GO terms that positively and negatively stimulate angiogenesis such as:
```
GO:0016525	angio	-	1	negative regulation of angiogenesis
GO:0045766	angio	+	1 	positive regulation of angiogenesis
GO:0043534	angio	+	1	blood vessel endothelial cell migration
GO:0043532	angio	-	1	angiostatin binding
```
With a defined SOI(s) and attributed GO terms, you can actually run the analysis and leave the other options at defaults. Other sections are explained in the following text.

#### Evidence code groups section
Evidence codes are three- or two-letter codes providing a specific level of proof for an annotation between a GO term and a specific gene. This section contains the whole hierarchy of possible evidence codes, grouped into several major evidence code groups (EGCs). This section only determines the possible EGCs and specific evidence codes, whereas the EGCs or specific evidence codes are selected in the _Settings_ section via the `evidence_codes` setting. 

Based on https://geneontology.org/docs/guide-go-evidence-codes/, there are the following 6 EGCs (noted with belonging evidence codes):
a. experimental evidence (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP) [experimental]
b. phylogenetically inferred evidence (IBA, IBD, IKR, IRD) [phylogenetic]
c. computational analysis evidence (ISS, ISO, ISA, ISM, IGC, RCA) [computational_analysis]
d. author statement evidence (TAS, NAS) [author_statement]
e. curator statement evidence (IC, ND) [curator_statement]
f. electronic annotation (IEA) [electronic]

Of important notice is that approximately 95% of Gene Ontology annotations are electronically inferred (IEA) and these are not checked by a human examiner.

This section exists to give user the option to add or exclude any evidence codes, should the GO evidence codes change in the future.
Each line contains two tab-separated elements:
- evidence code group name (e.g. author_statement)
- evidence codes (e.g. TAS,NAS) belonging to the group, along with their ECO identifiers (evidence code and identifier separated by underscore) as comma-separated values (e.g. TAS_ECO:0000304,NAS_ECO:0000303)

ECO evidence code identifiers can be found on https://wiki.geneontology.org/index.php/Guide_to_GO_Evidence_Codes and https://www.ebi.ac.uk/QuickGO/term/ECO:0000245.

WARNING: The evidence codes section MUST be specified before the settings section.

Example:
```
###evidence_code_groups
experimental	EXP_ECO:0000269,IDA_ECO:0000314,IPI_ECO:0000353,IMP_ECO:0000315,IGI_ECO:0000316,IEP_ECO:0000270,HTP_ECO:0006056,HDA_ECO:0007005,HMP_ECO:0007001,HGI_ECO:0007003,HEP_ECO:0007007
phylogenetic	IBA_ECO:0000318,IBD_ECO:0000319,IKR_ECO:0000320,IRD_ECO:0000321
computational_analysis	ISS_ECO:0000250,ISO_ECO:0000266,ISA_ECO:0000247,ISM_ECO:0000255,IGC_ECO:0000317,RCA_ECO:0000245
author_statement	TAS_ECO:0000304,NAS_ECO:0000303
curator_statement	IC_ECO:0000305,ND_ECO:0000307
electronic	IEA_ECO:0000501
```

#### Settings section
The settings section contains several settings, which are used to change the flow of the algorithm. 

**evidence_codes** is used to determine which annotations between GO terms and respective genes the algorithm will accept. GOReverseLookup will only accept genes annotated to input GO terms with any of the user-accepted evidence codes. 
- to accept all evidence codes belonging to a specific EGC, use a tilde operator in brackets `(~)`, e.g. `experimental(~)`
- to accept specific evidence codes belonging to an evidence group, specify them between the parentheses. If specific evidence codes are specified among parantheses, all non-specified evidence codes will be excluded. For example, to take into account only IC, but not ND, from curator_statement, use the following: `curator_statement(IC)`
- to exclude specific evidence codes, use an exclamation mark. All evidence not specified excluded evidence codes belonging to an EGC will still be included. To exclude only HEP and retain the rest of experimental evidence codes, use: `!experimental(HEP)`
- to merge multiple evidence code groups, supply them as comma-separated values. E.g.: `experimental(~),phylogenetic(~),computational_analysis(~),author_statement(TAS),curator_statement(IC),!electronic(~)`

Example evidence codes:
```
evidence_codes	experimental(~),phylogenetic(~),computational_analysis(~),author_statement(TAS),!curator_statement(ND),!electronic(~)
```

**pvalue** is the threshold _p_-value used to assess the statistical significance of a gene being involved in a target SOI. There are two possible cases of evaluation:

a) The user has defined an SOI and has attributed GO terms that both positively and negatively regulate the SOI. A gene is statistically significant if its p-value for the defined SOI stimulation/inhibition is less than the defined p-value threshold AND its p-value for the opposite SOI (inhibition/stimulation) is greater than the defined p-value threshold. It is advisable to also attribute GO terms that are opposite regulators of the defined target SOI in order to increase the credibility of the results.

b) The user has defined an SOI and has attributed GO terms only in one regulation direction (e.g. only stimulation or only inhibition). A gene is statistically significant if its p-value for the defined SOI is less than the defined p-value threshold.

**target_organism** is the target organism for which the statistical analysis is being performed. Organisms are represented with three identifiers (separated by vertical bars), which MUST be supplied for the program to correctly parse organism data: (1) organism label in lowercase (2) organism database and (3) organism NCBI taxon. For example, to select _Homo sapiens_ as the target organism, a researcher would specify:
```
target_organism	homo_sapiens|UniProtKB|NCBITaxon:9606
```

**ortholog_organisms** represent all homologous organisms, the genes of which are also taken into account during the scoring phase if they are found to have existing target organism orthologous genes. This feature has been enabled as a GO term can be associated with genes belonging to different organisms, which are indexed by various databases. The current model has been tested on the following orthologous organisms: _Rattus norvegicus_, _Mus musculus_, _Danio rerio_ and _Xenopus tropicalis_. Example:
```
ortholog_organisms	danio_rerio|ZFIN|NCBITaxon:7955,rattus_norvegicus|RGD|NCBITaxon:10116,mus_musculus|MGI|NCBITaxon:10090,xenopus_tropicalis|Xenbase|NCBITaxon:8364
```

**include_indirect_annotations**: if True, will increase the amount of annotations to a gene by the sum of all children GO terms of existing directly annotated GO terms to the gene. If False, will only count the direct annotations. This impacts the statistical relevance of genes during the scoring phase. Annotations from Gene Ontology between a GO term and a gene are directly annotated, but all children GO terms of the directly annotated term also infer the annotation. Consider the following tree:
```
GO:1901342 regulation of vasculature development
    - GO:0045765 regulation of angiogenesis
        - GO:0045766 positive regulation of angiogenesis <- gene Hipk2
            - GO:1905555 positive regulation of blood vessel branching
            - GO:1903672 positive regulation of sprouting angiogenesis
            - GO:0035470 positive regulation of vascular wound healing
```
Gene Hipk2 is directly annotated to GO:0045766. The children annotations also infer the annotation (GO:1905555, GO:1903672, GO:0035470), but not the parent annotation (GO:1901342). 

**goterm_gene_query_timeout** is the timeout it takes when querying genes annotated to GO terms. If specifying very vague GO terms (such as `regulation of gene expression`, which has ~25 million annotations, a query might fail due to a request taking too long to complete or, which is a more severe error due to its covertness, a query might return an incomplete list of genes associated with a GO term. As a rule of thumb, we discourage the usage of such vague GO terms. A default 240-second timeout ensures that all GO terms approximately with a few million annotations are fetched correctly from the GO servers.

**goterm_gene_query_max_retries** is the maximum number of retries sent to the GO servers before dropping a GO term and assigning it with an empty list of associated genes.

**gorth_ortholog_refetch**
We implemented a gOrth batch ortholog query (https://biit.cs.ut.ee/gprofiler/orth), which speeds up the total runtime of the program. The function attempts to find orthologs to genes in a single batch request. If 'gorth_ortholog_refetch' is True, then the genes for which orthologs were not found will be re-fetched using alternative Ensembl calls. If 'gorth_ortholog_refetch' is False, then the genes for which orthologs were not found will not be queried for orthologs again.

**gorth_ortholog_fetch_for_indefinitive_orthologs**
The gOrth batch query implementation can return the following options:
- multiple orthologous genes (these are called "indefinitive orthologs")
- a single orthologous gene (called a "definitive ortholog")
- no orthologous genes.

In our asynchronous Ensembl ortholog query pipeline implementation, when multiple orthologous genes are returned from Ensembl, the orthologous gene with the highest percentage identity (percentage identity of amino-acid sequence between the gene and the target organism orthologous gene) is selected as the best ortholog and is assigned as the true ortholog to the input gene. However, gOrth has currently (10_29_2023) no option to return the "best" orthologous gene, neither it has the option to exclude obsolete ortholog gene ids (confirmed by the gProfiler team via an email conversation). Therefore, it is advisable to keep the gorth_ortholog_fetch_for_indefinitive_orthologs to True, so that indefinitive orthologs are discarded from the gOrth ortholog query and are instead fetched by the asynchronos pipeline, which can select the best ortholog for the input gene. Having this setting set to False will choose, in the case of indefinitive orthologs, the first returned ortholog id from the gOrth query, but with no guarantees that this ortholog id is not obsolete.

**fisher_test_use_online_query**
It is highly advisable to leave this setting set to False, otherwise, the timing of the scoring phase might severely be extended (into days, if not weeks).

**uniprotkb_genename_online_query**: When querying all genes associated to a GO Term, Gene Ontology returns UniProtKB identified genes (amongst others, such as ZFIN, Xenbase, MGI, RGD). During the algorithm, gene name has to be determined. It can be obtained via two pathways:

- online pathway, using UniProtAPI
- offline pathway, using the GO Annotations File

During testing, it has been observed that the offline pathway usually results in more gene names found, besides being much faster. Thus, it is advisable to leave this setting set to False, both to increase speed and accuracy. If it is set to True, then gene names will be queried from the UniProtKB servers.


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

# GOReverseLookup

[![PyPI package](https://img.shields.io/badge/pip%20install-goreverselookup-brightgreen)](https://pypi.org/project/goreverselookup/) [![version number](https://img.shields.io/github/v/release/MediWizards/GOReverseLookup)](https://github.com/MediWizards/GOReverseLookup/releases) [![Actions Status](https://img.shields.io/github/actions/workflow/status/MediWizards/GOReverseLookup/test_on_push.yml)](https://github.com/MediWizards/GOReverseLookup/actions/workflows/test_on_push.yml) [![License](https://img.shields.io/github/license/MediWizards/GOReverseLookup)](https://github.com/MediWizards/GOReverseLookup/blob/main/LICENSE)

**GOReverseLookup** is a Python package designed for Gene Ontology Reverse Lookup. It serves the purpose of identifying statistically significant genes within a set or a cross-section of selected Gene Ontology Terms. Researchers need only define their own states of interest (SOIs), and select GO terms must be attributed as either positive or negative regulators of the chosen SOIs. For more information regarding the creation of the input file for the program, refer to the _Input file_ section. Once the input file is completed, the GOReverseLookup program can be started. Once the algorithm is completed, the program saves statistically significant genes in a standalone file.

For example, if researchers were interested in the _angiogenesis_ SOI, then an attributed group of GO terms as positive regulators of _angiogenesis_ might have been defined using the following GO terms:
- GO:1903672 positive regulation of sprouting angiogenesis
- GO:0001570 vasculogenesis
- GO:0035476 angioblast cell migration

And negative regulators of the _angiogenesis_ SOI might have been defined as the following group:
- GO:1903671 negative regulation of sprouting angiogenesis
- GO:1905554 negative regulation of vessel branching
- GO:0043537 negative regulation of blood vessel endothelial cell migration

If a researcher defines the _target process_ as positive regulation of a desired SOI (in our case _angiogenesis_), then GOReverseLookup finds all genes statistically relevant for the group of GO terms defined as positive regulators of _angiogenesis_ (p < 0.05), while excluding any genes determined to be statistically significant (_p_ < 0.05) in the opposing process (in our case, negative regulation of _angiogenesis_). _P_-value threshold can also be manually set by the user.


## Getting Started

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

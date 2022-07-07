# greatpy

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/workflow/status/ilibarra/greatpy/Test/main
[link-tests]: https://github.com/theislab/greatpy/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/greatpy

Implementation of GREAT in Python

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

### What is great : 
Great is a bioinformatics tool to analyze cis-regulatory regions of the genome. This method can be used to assign probable biological meanings to unannotated genes based on the annotations carried by neighboring genes. This tool allows going from peaks in .bed format to enrichment ontology terms associated with these peaks using several statistical tests. 
This package is strongly inspired by [great][great_article] allowing Helmholtz to have a stable and perennial version of the package.

<img align="center" src="./img/great_tool.png?raw=true">


### What can you do with greatpy : 
* Translate a genetic file in .bed format and containing the following information: chromosome number, start position on the chromosome, end position, gene name and tss. Into a regulatory region file that can then be used in the great : 
```python 
regdom = great.tl.create_regdom(
    tss_file=path_of_the_file,
    chr_sizes_file="path_of_the_file",
    association_rule="Basalplusextention",
    out_path=path_output
    )
```
The association's rules could be : 
* Basalplusextention 
* OneCloset 
* TwoCloset

Documentation available [here][association_rules]: 

* Analyzes the significance of proximal and distal cis-regulatory regions in the genome. To do this: 
```python 
res = great.tl.enrichment(
    test=path_of_genomic_region_to_test,
    regdom_file=path_of_regdom_file,
    chr_size_file=path_each_chromosome_size,
    annotation=path_of_the_csv_file_of_ontologies,
```
Several arguments can be added to this function such as : 
* binom (default true): should the binomial p-value be calculated 
* hypergeom (default true): should the hypergeometric p-value be computed 
* alpha (default 0.05) : significance threshold 
* correction (default ("fdr",0.05)): Which correction should be applied, correction[0] = "bonferroni" | "fdr" | 0 ; correction[1] = initial significance threshold. 
* sort_by: by which column the results should be sorted 

Several examples of uses can be found in the notebook part of the package: 
* For the create_regdom option: [explanatory notebook][notebook1]
* For the enrichment function: [explanatory notebook] []

## Installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`\_.

There are several alternative options to install greatpy:

<!--
1) Install the latest release of `greatpy` from `PyPI <https://pypi.org/project/greatpy/>`_:

```bash
pip install greatpy
```
-->

1. Install the latest development version:

```bash
pip install git+https://github.com/ilibarra/greatpy.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

> t.b.a

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/ilibarra/greatpy/issues
[changelog]: https://greatpy.readthedocs.io/latest/changelog.html
[link-docs]: https://greatpy.readthedocs.io
[link-api]: https://greatpy.readthedocs.io/latest/api.html
[great_article]: https://www.nature.com/articles/nbt.1630
[association_rules]: https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules
[notebook1]: https://github.com/theislab/greatpy/blob/main/notebooks/01_create_regdom.ipynb
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

### What is GREAT :

GREAT (Genomic Regions Enrichment of Annotations Tool) is a bioinformatics tool, this method enables to associate genetic regions to the most probable GO terms.

### How to use greatpy :

This package is strongly inspired by [GREAT][great_article] allowing Helmholtz to have a stable, perennial and updated version of the package.

[GREAT figure][great_figure] issue from [GREAT article][great_article]

#### <ins>1. Create regulatory domain from tss</ins>

-   Translate a genetic file in .bed format and containing the following information: chromosome number, start position on the chromosome, end position, gene name and TSS.
    Into a regulatory region file that can then be used in greatpy:

```python
regdom = greatpy.tl.REGDOM.create_regdom(
    tss_file=path_of_the_file,
    chr_sizes_file="path_of_the_file",
    association_rule="Basalplusextention",
    out_path=path_save_output,
)
```

The [association rules][association_rules] parameters could be :

-   `Basalplusextention`
-   `OneCloset`
-   `TwoCloset`

<p align="center">
  <img src="./sketch/association_rule.jpg?raw=true" style="width:75%">
</p>

#### <ins>2. Get enrichment of GO term in the tests genomics regions</ins>

-   Analyzes the significance of proximal and distal cis-regulatory regions in the genome. To do this:

```python
res = greatpy.tl.GREAT.enrichment(
    test=path_of_genomic_region_to_test,
    regdom_file=path_of_regdom_file,
    chr_size_file=path_each_chromosome_size,
    annotation=path_of_the_csv_file_of_ontologies,
)
```

Several arguments can be added to this function such as :

-   `binom` (default True): should the binomial p-value be calculated?
-   `hypergeom` (default True): should the hypergeometric p-value be computed?

It is then possible to apply a Bonferroni and/or FDR correction to the found p-values:

```python
res = greatpy.tl.GREAT.enrichment(
    test_file=path_or_dataframe_of_genomic_region_to_test,
    regdom_file=path_or_dataframe_of_regdom_file,
    chr_size_file=path_or_dataframe_each_chromosome_size,
    annotation_file=path_or_dataframe_of_the_csv_file_of_ontologies,
)
great.tl.GREAT.set_fdr(res, alpha=0.05)
great.tl.GREAT.set_bonferroni(res, alpha=0.05)
```

#### <ins>3. Plot</ins>

It is also possible to create several types of plots :

-   Number of genetic associations per genomic region
-   Distance to the associated gene TSS for each genomic region studied
-   Absolute distance to the associated gene TSS for each genomic region studied

```python
fig, ax = plt.subplots(1, 3, figsize=(30, 8))
greatpy.pl.graph_nb_asso_per_peaks(
    path_or_dataframe_of_genomic_region_to_test, path_or_dataframe_of_regdom_file, ax[0]
)
greatpy.pl.graph_dist_tss(
    path_or_dataframe_of_genomic_region_to_test, path_or_dataframe_of_regdom_file, ax[1]
)
greatpy.pl.graph_absolute_dist_tss(
    path_or_dataframe_of_genomic_region_to_test, path_or_dataframe_of_regdom_file, ax[2]
)
plt.show()
```

<p align="center">
  <img src="./sketch/plot1.png?raw=true">
</p>

-   Dotplot of the enrichment GO term in the genomic test region

```python
plot = enrichment_df.rename(columns={"binom_p_value": "p_value", "go_term": "name"})
plt.figure(figsize=(10, 10))
great.pl.plot_enrich(plot)
```

<p align="center">
  <img src="./sketch/dotplot.png?raw=true" style="width:75%">
</p>

-   Dotplot of the enrichment GO terms in multiple genomic test regions

```python
test = ["name_bindome_biosample_1", "name_bindome_biosample_2", "..."]
tmp_df = great.tl.GREAT.enrichment_multiple(
    tests=test,
    regdom_file="../data/human/hg38/regulatory_domain.bed",
    chr_size_file="../data/human/hg38/chr_size.bed",
    annotation_file="../data/human/ontologies.csv",
    binom=True,
    hypergeom=True,
)
```

<p align="center">
  <img src="./sketch/multidot.png?raw=true" alt="dotplot of multi sample" width="300" height="400">
</p>

Several examples of uses can be found in the notebook part of the package:

-   For the create_regdom option: [notebook][notebook1]
-   For the enrichment function: [notebook][notebook2]
-   Plot : [notebook][notebook3]

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
pip install git+https://github.com/theislab/greatpy.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

For cite greatpy :

```bibtex
@software{greatpy,
  author = {Ibarra, Mauger-Birocheau},
  doi = {},
  month = {},
  title = {{greatpy}},
  url = {https://github.com/theislab/greatpy},
  year = {2022}
}
```

## References

```bibtex
@article{GREAT,
  author   = {McLean, C.
              and Bristor, D.
              and Hiller, M. et al.},
  title    = {GREAT improves functional interpretation of cis-regulatory regions},
  journal  = {Nat Biotechnol},
  year     = {2010},
  month    = {May},
  day      = {02},
  volume   = {28},
  number   = {495},
  pages    = {501},
  doi      = {10.1038/nbt.1630},
  url      = {https://doi.org/10.1038/nbt.1630}
}
```

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/ilibarra/greatpy/issues
[changelog]: https://greatpy.readthedocs.io/latest/changelog.html
[link-docs]: https://greatpy.readthedocs.io
[link-api]: https://greatpy.readthedocs.io/latest/api.html
[great_article]: https://www.nature.com/articles/nbt.1630
[great_figure]: https://www.nature.com/articles/nbt.1630/figures/1
[association_rules]: https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules
[notebook1]: https://github.com/theislab/greatpy/blob/main/notebooks/01_create_regdom.ipynb
[notebook2]: https://github.com/theislab/greatpy/blob/main/notebooks/02_binom_vs_hypergeom.ipynb
[notebook3]: https://github.com/theislab/greatpy/blob/main/notebooks/07_plot.ipynb

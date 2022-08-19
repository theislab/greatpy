# greatpy

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/workflow/status/ilibarra/greatpy/Test/main
[link-tests]: https://github.com/theislab/greatpy/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/greatpy

Implementation of GREAT in Python

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

## Notebook

|   Information             |   link                    |
| ------------------------- | ------------------------- |
|   Create regdom           | [notebook][notebook1]     |
|   enrichment              | [notebook][notebook2]     |
|   plot                    | [notebook][notebook3]     |
|   Comparaison with GREAT  | [notebook][notebook4]     |

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

### What is GREAT :

GREAT (Genomic Regions Enrichment of Annotations Tool) is a bioinformatics tool, this method enables to associate genetic regions to the most probable GO terms.

### How to use greatpy :

This package is strongly inspired by [GREAT][great_article] allowing Helmholtz to have a stable, perennial and updated version of the package.

[GREAT figure][great_figure] issue from [GREAT article][great_article]

#### <ins>1. Create regulatory domain from tss</ins>

-   Translate a genetic file in `.bed` format and containing the following information: 
  -   TSS file should have the following columns :`\t` `chromosome_number` `\t` `position` `\t` `strand` `\t` `gene_name`.
  -   Chromosome size file should have the following columns :`\t` `chromosome_number` `\t` `chromosome_size`.

```python
regdom = greatpy.tl.REGDOM.create_regdom(
    tss_file=Input_TSS_path, # eg : "../data/human/hg38/tss.bed"
    chr_sizes_file=Input_chromosome_size_path, # eg : "../data/human/hg38/chr_size.bed"
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

-   Analyzes the significance of proximal and distal cis-regulatory regions in the genome.
-   Some files should be used as input : 
    -   test file should have the following columns :`\t` `chr` `\t` `chr_start` `\t` `chr_end`.
    -   regulatory domain file should have the following columns :`chr` `\t` `chr_start` `\t` `chr_end` `\t` `name` `\t` `tss	strand`
    -   chromosome size file should have the following columns :`\t` `chromosome_number` `\t` `chromosome_size`.
    -   annotation file should have the following columns :`\t` `ensembl` `\t` `id` `\t` `name` `\t` `ontology.group` `\t` `gene.name` `\t` `symbol`

```python
res = greatpy.tl.GREAT.enrichment(
    test_file=Input_path_or_df, # eg : "../data/tests/test_data/input/10_MAX.bed"
    regdom_file=regdom_path_or_df, # eg : "../data/human/hg38/regdom.bed"
    chr_size_file=chromosome_size_path_or_df, # eg : "../data/human/hg38/chr_size.bed"
    annotation_file=annotation_path_or_df, # eg : "../data/human/ontologies.csv"
)
```

Several arguments can be added to this function such as :

-   `binom` (default True): should the binomial p-value be calculated?
-   `hypergeom` (default True): should the hypergeometric p-value be computed?

It is then possible to apply a Bonferroni and/or FDR correction to the found p-values:

```python
great.tl.GREAT.set_fdr(res, alpha=0.05)
great.tl.GREAT.set_bonferroni(res, alpha=0.05)
```

#### <ins>3. Plot</ins>

##### 1 genomic distribution of data

-   Number of genetic associations per genomic region
-   Distance to the associated gene TSS for each genomic region studied
-   Absolute distance to the associated gene TSS for each genomic region studied

```python
fig, ax = plt.subplots(1, 3, figsize=(30, 8))
greatpy.pl.graph_nb_asso_per_peaks(
    Input_path_or_df, # eg : "../data/tests/test_data/input/10_MAX.bed"
    regdom_path_or_df, # eg : "../data/human/hg38/regdom.bed"
    ax[0]
)
greatpy.pl.graph_dist_tss(
    Input_path_or_df, # eg : "../data/tests/test_data/input/10_MAX.bed"
    regdom_path_or_df, # eg : "../data/human/hg38/regdom.bed"
    ax[0]
)
greatpy.pl.graph_absolute_dist_tss(
    Input_path_or_df, # eg : "../data/tests/test_data/input/10_MAX.bed"
    regdom_path_or_df, # eg : "../data/human/hg38/regdom.bed"
    ax[0]
)
plt.show()
```

<p align="center">
  <img src="./sketch/plot1.png?raw=true">
</p>

##### 2 Enrichments by GO terms (dotplot) - one input

```python
plot = enrichment_df.rename(columns={"binom_p_value": "p_value", "go_term": "name"})
plt.figure(figsize=(10, 10))
great.pl.plot_enrich(plot)
```

<p align="center">
  <img src="./sketch/dotplot.png?raw=true" style="width:75%">
</p>

#### 3 Enrichments by GO terms (dotplot) - multiple inputs

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

## Note

Both types of tests (binomial and hypergeometric) performed may be susceptible to certain biases of which one must be aware to analyze the results with a critical mind.

-   The hypergeometric test may be biased by the size of the regulatory domains of the genes since isolated genes have very large regulatory domains and are therefore more likely to generate false positives.
-   The binomial test can also be biased if a large number of genomic regions to be tested are associated with a small set of genes that can also generate false positives.

But these biases are partially compensated between each of the tests the binomial test reduces the hypergeometric bias by taking into account exactly the size of the regulatory domains of the genes and the hypergeometric test compensates for the bias of the binomial test by counting each gene only once.
The two types of tests are complementary and must be analyzed together to conclude.

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

```bibtex
@Manual{rGREAT,
  title = {rGREAT: GREAT Analysis - Functional Enrichment on Genomic Regions},
  author = {Zuguang Gu},
  year = {2022},
  note = {https://github.com/jokergoo/rGREAT, http://great.stanford.edu/public/html/},
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
[notebook4]: https://greatpy.readthedocs.io/en/latest/notebooks/03_great_vs_greatpy.html

# Example

## Code snippet

### Create regulatory domain

greatpy can allow you to create regulatory domains with a TSS.bed and chromosome_size.bed files with the function `greatpy.tl.create_regdom()`. You can show an example of this function with the following code:

```
import greatpy as great
great.tl.REGDOM.create_regdom(
    tss_file = '/path/to/tss.bed',
    chrom_size_file = '/path/to/chrom_size.bed',
    association_rule = "Basalplusextention_or_OneCloset_ore_TwoCloset",
    out_path = "None_or_path/to/save/output"
)
```

### Get enrichment

greatpy can allow you to compute the GO term enrichment on a set of the chromosomic region on a bed format with `greatpy.tl.enrichment()`. You can show an example of this function with the following code:

```
import greatpy as great
enrichment = great.tl.enrichment(
    test_file = df_or_path_to_test_file,
    regdom_file = df_or_path_to_regdomfile_file, # could be create with great.tl.REGDOM.create_regdom()
    chr_size_file = df_or_path_to_size_file,
    annotation_file = df_or_path_to_ontology_annotation_file,
    binom = True,
    hypergeom = True,
    )
```

After the calculation, it is possible to apply corrections to the p_values, two methods of corrections are possible:

-   Bonferroni correction :

```
enrichment = great.tl.set_fdr(enrichment)
```

-   FDR correction :

```
enrichment = great.tl.set_bonferroni(enrichment)
```

It is also possible to apply a threshold on one of the columns to reduce the table

```
enrichment = great.tl.set_threshold(enrichment,colname = "column_to_apply_the_threshold",alpha = 0.05)
```

### plot the results

Several types of plots can be made:

```
import greatpy as great
import matplotlib.pyplot as plt
```

#### Scatter plot

```
great.pl.scatterplot(
    enrichment,
    x = "colname",
    y = "colname",
    title = "title of the plot",
    xlabel = "x_label",
    ylabel = "y label"
)
plt.show()
```

```{image} _static/output_images/scatterplot.png
:width: 500px
```

#### Graph of the number of associations per peak

```
fig,ax = plt.subplots(1,3,figsize = (30,8))
great.pl.graph_nb_asso_per_peaks(test,regdom,ax[0])
great.pl.graph_dist_tss(test,regdom,ax[1])
great.pl.graph_absolute_dist_tss(test,regdom,ax[2])
plt.show()
```

```{image} _static/output_images/plot1.png

```

#### Dotplot showing the enrichment of the GO terms

```
plot = enrichment.rename(columns = {"binom_p_value": "p_value", "go_term": "name"})
great.pl.plot_enrich(plot)
```

```{image} _static/output_images/dotplot.png
:width: 600px
```

```
test = [
    "SRF:Ishikawa,A-673-clone-Asp114,K-562,MCF-7,Hep-G2",
    "MAX:K-562,WA01,HeLa-S3", "BACH1:A-549,GM12878",
    "CDK9:A-375,MM1-S,MV4-11,P493-6,BT-474,HEK293T",
    "GATA1:erythroblast,HUDEP-2,K-562",
    "IKZF1:K-562,GM12878,HSPC",
    "SP1:liver,A-375,Hep-G2,HEK293,GM12878,A-549,K-562,HEK293T,WA01",
    "TCF7:Hep-G2,GM12878,K-562",
    "ZBTB40:MCF-7,Hep-G2,GM12878",
    "AFF1:MV4-11,K-562"
    ]

results = great.tl.enrichment_multiple(
    tests = test,
    regdom_file = "../data/human/hg38/regulatory_domain.bed",
    chr_size_file = "../data/human/hg38/chr_size.bed",
    annotation_file = "../data/human/ontologies.csv",
    binom = True,
    hypergeom = True,
    )

fig = plt.figure(figsize = (15, 12))
p_val,odd_ratio,df = great.pl.dotplot_multi_sample(
    results,
    fig = fig,
    show_term_name = True,
    term_name_nchars = 20
    )
```

```{image} _static/output_images/multidot.png
:width: 600px
```

## Notebook example

```{toctree}
:hidden: true
:maxdepth: 1

notebooks/01_create_regdom.ipynb
notebooks/02_binom_vs_hypergeom.ipynb
notebooks/03_great_vs_greatpy.ipynb
notebooks/04_hypergeom_vs_gprofiler.ipynb
notebooks/05_rgreat_online_vs_local_vs_greatpy.ipynb
notebooks/06_peaks_using_bindome.ipynb
notebooks/07_plot.ipynb
```

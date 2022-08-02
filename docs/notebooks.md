# Example

## Code snippet

### Create regulatory domain

greatpy can allow you to create regulatory domains with a TSS.bed and chromosome_size.bed files with the function `greatpy.tl.REGDOM.create_regdom()`. You can show an example of this function with the following code:

```
import greatpy as gp
gp.tl.REGDOM.create_regdom(
    tss_file='/path/to/tss.bed',
    chrom_size_file='/path/to/chrom_size.bed',
    association_rule="Basalplusextention_or_OneCloset_ore_TwoCloset",
    out_path=None_or_path/to/save/output
)
```

### Get enrichment

greatpy can allow you to compute the GO term enrichment on a set of the chromosomic region on a bed format with `greatpy.tl.GREAT.enrichment()`. You can show an example of this function with the following code:

```
import greatpy as gp
enrichment = gp.tl.GREAT.enrichment(
    test_file=df_or_path_to_test_file,
    regdom_file=df_or_path_to_regdomfile_file, # could be create with gp.tl.REGDOM.create_regdom()
    chr_size_file=df_or_path_to_size_file,
    annotation_file=df_or_path_to_ontology_annotation_file,
    binom=True,
    hypergeom=True,
    )
```

After the calculation, it is possible to apply corrections to the p_values, two methods of corrections are possible:

-   Bonferroni correction :

```
gp.tl.GREAT.set_fdr(enrichment)
```

-   FDR correction :

```
gp.tl.GREAT.set_bonferroni(enrichment)
```

It is also possible to apply a threshold on one of the columns to reduce the table

```
gp.tl.GREAT.set_threshold(enrichment,colname="column_to_apply_the_threshold",alpha=0.05)
```

### plot the results

Several types of plots can be made:

```
import greatpy as gp
import matplotlib.pyplot as plt
```

#### Scatter plot

```
gp.pl.scatterplot(
    enrichment,
    x="colname",
    y="colname",
    title="title of the plot",
    xlabel="x_label",
    y label="y label"
)
plt.show()
```

```{image} _static/output_images/scatterplot.png
:width: 200px
```

#### Graph of the number of associations per peak

```
fig,ax = plt.subplots(1,3,figsize=(30,8))
gp.pl.graph_nb_asso_per_peaks(test,regdom,ax[0])
gp.pl.graph_dist_tss(test,regdom,ax[1])
gp.pl.graph_absolute_dist_tss(test,regdom,ax[2])
plt.show()
```

```{image} _static/output_images/plot1.png

```

#### Dotplot showing the enrichment of the GO terms

```
plot = enrichment.rename(columns={"binom_p_value" : "p_value", "go_term":"name"})
gp.pl.plot_enrich(plot)
```

```{image} _static/output_images/dotplot.png
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
notebooks/06_peaks_using_binodome.ipynb
notebooks/07_plot.ipynb
```

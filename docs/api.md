```{eval-rst}
.. module:: greatpy
```

```{eval-rst}
.. automodule:: greatpy
   :noindex:
```

# API

```
import greatpy as gp
```

## Tools : `tl`

```{eval-rst}
.. module:: greatpy.tl
```

```{eval-rst}
.. currentmodule:: greatpy
```

### Create regulatory domain

```{eval-rst}
.. autosummary::
    :toctree: generated/

    tl.regdom.Regdom.create_regdom
```

### greatpy computation

#### Main functions

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.great.Great.enrichment
    tl.great.Great.enrichment_multiple
```

#### Additional functions

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.great.Great.loader
    tl.great.Great.set_bonferroni
    tl.great.Great.set_fdr
    tl.great.Great.set_threshold
```

### utils additional tools

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.utils.get_nb_asso_per_region
    tl.utils.get_dist_to_tss
    tl.utils.online_vs_local_vs_greatpy_comparison
    tl.get_association
    tl.len_regdom
    tl.number_of_hits
    tl.get_binom_pval
    tl.hypergeom_pmf
    tl.hypergeom_cdf
```

## Plotting : `pl`

```{eval-rst}
.. module:: greatpy.pl
.. currentmodule:: greatpy

.. autosummary::
    :toctree: generated/

    pl.scatterplot
    pl.graph_nb_asso_per_peaks
    pl.graph_dist_tss
    pl.graph_absolute_dist_tss
    pl.plot_enrich
    pl.make_bubble_heatmap
    pl.dotplot_multi_sample
    pl.get_all_comparison
```

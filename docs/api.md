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

    tl.REGDOM.create_regdom
```

### greatpy computation
```{eval-rst}

.. autosummary::
    :toctree: generated

    tl.GREAT.loader
    tl.GREAT.enrichment
    tl.GREAT.set_bonferroni
    tl.GREAT.set_fdr
    tl.GREAT.set_threshold
```
### utils side tools
```{eval-rst}
.. currentmodule:: greatpy

.. autosummary::
    :toctree: generated

    tl.utils.get_nb_asso_per_region 
    tl.utils.get_dist_to_tss
    tl.get_association
    tl.len_regdom 
    tl.number_of_hit
    tl.betacf
    tl.betai
    tl.get_binom_pval 
    tl.hypergeom_pmf 
    tl.hypergeom_cdf
```

## Plotting : `tl`
```{eval-rst}
.. currentmodule:: greatpy

.. autosummary::
    :toctree: generated

    pl.scatterplot
    pl.graph_nb_asso_per_peaks
    pl.graph_dist_tss
    pl.graph_absolute_dist_tss
    pl.plot_enrich
```

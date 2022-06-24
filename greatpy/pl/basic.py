from anndata import AnnData
import pandas as pd 
from numpy import log
from seaborn import scatterplot

def basic_plot(adata: AnnData) -> int:
    """Generate a basic plot for an AnnData object."""
    print("Import matplotlib and implement a plotting function here.")
    return 0

def scatterplot_of_p_val(Great_df:pd.DataFrame,colname_x,colname_y):
    Great_df = Great_df.dropna()
    Great_df = Great_df.astype({"Binom p-value":"float64","Hypergeom p-value":"float64"})
    Great_df[f"-log({colname_x})"] = -log(Great_df[colname_x])
    Great_df[f"-log({colname_y})"] = -log(Great_df[colname_y])

    scatterplot(data=Great_df,x=f"-log({colname_x})",y=f"-log({colname_y})")
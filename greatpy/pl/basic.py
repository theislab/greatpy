from anndata import AnnData
import pandas as pd 
from numpy import log
from seaborn import scatterplot as sp

def basic_plot(adata: AnnData) -> int:
    """Generate a basic plot for an AnnData object."""
    print("Import matplotlib and implement a plotting function here.")
    return 0

def scatterplot(great_df:pd.DataFrame,colname_x,colname_y,minus_log10=True):
    great_df = great_df.dropna()
    great_df = great_df.astype({colname_x:"float64",colname_y:"float64"})
    if minus_log10 :
        great_df[f"-log({colname_x})"] = -log(great_df[colname_x])
        great_df[f"-log({colname_y})"] = -log(great_df[colname_y])
        sp(data=great_df,x=f"-log({colname_x})",y=f"-log({colname_y})")
    else : 
        sp(data=great_df,x=colname_x,y=colname_y)
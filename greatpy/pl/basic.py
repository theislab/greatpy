from anndata import AnnData
import pandas as pd 
from numpy import log
from seaborn import scatterplot as sp,barplot as bar 
import matplotlib.pyplot as plt

def basic_plot(adata: AnnData) -> int:
    """Generate a basic plot for an AnnData object."""
    print("Import matplotlib and implement a plotting function here.")
    return 0

def scatterplot(great_df:pd.DataFrame,colname_x,colname_y,title:str="",minus_log10=True,ax=None):
    great_df = great_df.dropna()
    great_df = great_df.astype({colname_x:"float64",colname_y:"float64"})
    if minus_log10 :
        great_df[f"-log({colname_x})"] = -log(great_df[colname_x])
        great_df[f"-log({colname_y})"] = -log(great_df[colname_y])
        sp(data=great_df,x = f"-log({colname_x})",y = f"-log({colname_y})",ax = ax).set_title(title)
    else : 
        sp(data=great_df,x=colname_x,y=colname_y,ax=ax).set_title(title)

###################################################################
########################### unfinish ##############################
###################################################################
def graph_nb_asso_per_peaks(dict_nb_asso_per_peaks:dict) :
    pass

def graph_dist_tss(dist_tss:dict, abs:bool = False) : 
    """Plot the distribution of TSS distance."""
    if not abs : 
        res = {"<-500":[0],"-500:-50":[0],"-50:-5":[0],"-5:0":[0],"0:5":[0],"5:50":[0],"50:500":[0],">500":[0]}
        for dist in dist_tss.values() :
            if dist < -500000 :
                res["<-500"][0] += 1
            elif dist < -50000 :
                res["-500:-50"][0] += 1
            elif dist < -5000 :
                res["-50:-5"][0] += 1
            elif dist < 0 :
                res["-5:0"][0] += 1
            elif dist < 5000 :
                res["0:5"][0] += 1
            elif dist < 50000 :
                res["5:50"][0] += 1
            elif dist < 500000 :
                res["50:500"][0] += 1
            else :
                res[">500"][0] += 1
        res = pd.DataFrame(res,index=["count"]).transpose().reset_index(drop=False).rename(columns={"index":"dist"})
        bar(x = "dist", y = "count",data = res).title("Distribution of distance to TSS")
        plt.show()
    else :
        res = {"0:5":[0],"5:50":[0],"50:500":[0],">500":[0]}
        for dist in dist_tss.values() :
            dist = abs(dist)
            if dist < 5000 :
                res["0:5"][0] += 1
            elif dist < 50000 :
                res["5:50"][0] += 1
            elif dist < 500000 :
                res["50:500"][0] += 1
            else :
                res[">500"][0] += 1
        res = pd.DataFrame(res,index=["count"]).transpose().reset_index(drop=False).rename(columns={"index":"dist"})
        bar(x = "dist", y = "count",data = res).title("Distribution of absolute distance to TSS")
        plt.show()
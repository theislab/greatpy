import pandas as pd 
from numpy import log
from seaborn import scatterplot as sp,barplot as bar 
import greatpy as gp

def scatterplot(great_df:pd.DataFrame,colname_x:str,colname_y:str,title:str="",minus_log10=True,ax=None) -> None :
    """
    This function is used to create a scatterplot from a 
    pandas dataframe between two columns. 
    A logarithmic scale can be used. 

    Parameters
    ----------
    great_df : pd.DataFrame
        Output of the greatpy.tl.GREAT.enrichment function
    colname_x : str 
        Name of the column to be used as x axis
    colname_y : str
        Name of the column to be used as y axis
    title : str
        Title of the plot
    minus_log10 : bool
        If True, the logarithmic scale is used
    ax : 
        Define the position of the plot in a figure 

    Returns
    -------
    None

    Exemples 
    --------
    Example available here: https://github.com/theislab/greatpy/blob/main/notebooks/02_binom_vs_hypergeom.ipynb

    """
    great_df = great_df.dropna()
    great_df = great_df.astype({colname_x:"float64",colname_y:"float64"})
    if minus_log10 :
        great_df[f"-log({colname_x})"] = -log(great_df[colname_x])
        great_df[f"-log({colname_y})"] = -log(great_df[colname_y])
        sp(data=great_df,x = f"-log({colname_x})",y = f"-log({colname_y})",ax = ax).set_title(title)
    else : 
        sp(data=great_df,x=colname_x,y=colname_y,ax=ax).set_title(title)


def graph_nb_asso_per_peaks(test:str or pd.DataFrame,regdom:str or pd.DataFrame,ax=None,color=None) -> None :
    """
    This function creates a barplot representing the 
    percentage of peaks for all possible association numbers  

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame 
        Regulatory domain of all genes in the genome 
    ax : 
        Define the position of the plot in a figure 

    Returns
    -------
    None

    Exemples 
    --------
    Example available here: https://github.com/theislab/greatpy/blob/main/notebooks/07_plot.ipynb   
    """
    nb_asso_per_peaks = gp.tl.get_nb_asso_per_region(test,regdom)

    nb = {
        "number" : [],
        "number_genes" : [],
        "percentage" : [],
    }
    for i in list(set(nb_asso_per_peaks.values())) :
        nb["number"].append(i)
        nb["number_genes"].append(list(nb_asso_per_peaks.values()).count(i))
        nb["percentage"].append(round(( list(nb_asso_per_peaks.values()).count(i) / len(nb_asso_per_peaks.keys()) ) * 100))
    nb = pd.DataFrame(nb,columns=["number","number_genes","percentage"],index=nb["number"])

    g = bar(data = nb,x="number",y="percentage",ax=ax,color=color)
    g.set_title("Number of associated genes per region")
    g.set_xlabel("Number of associated genes per region")
    g.set_ylabel("Genomic region (%)")

    for i in range(nb.shape[0]):  
        x = nb.iloc[i]["number"]
        y = nb.iloc[i]["percentage"]
        g.text(x = x - 0.06,y=y + 1,s = nb.number_genes[0])

def graph_dist_tss(test:str or pd.DataFrame,regdom:str or pd.DataFrame,ax=None,color="#325fa8") -> None : 
    """
    This function allows the creation of a barplot of the distance 
    between the peaks and the TSS of the associated gene(s). 

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame 
        Regulatory domain of all genes in the genome 
    ax : 
        Define the position of the plot in a figure

    Returns
    -------
    None

    Exemples 
    --------
    Example available here: https://github.com/theislab/greatpy/blob/main/notebooks/07_plot.ipynb 
   
    """
    res = {"<-500": [0],"-500:-50": [0],"-50:-5": [0],"-5:0": [0],"0:5": [0],"5:50": [0],"50:500": [0],">500": [0]}
    nb = 0

    dist = gp.tl.get_dist_to_tss(test,regdom)
    for i in dist.values() : 
        for j in i : 
            if j < -500000 : 
                res["<-500"][0] += 1
            elif j < -50000 : 
                res["-500:-50"][0] += 1
            elif j < -5000 :
                res["-50:-5"][0] += 1
            elif j < 0 :
                res["-5:0"][0] += 1
            elif j < 5000 :
                res["0:5"][0] += 1
            elif j < 50000 :
                res["5:50"][0] += 1
            elif j < 500000 :
                res["50:500"][0] += 1
            else :
                res[">500"][0] += 1
            nb += 1
    df = pd.DataFrame(res).transpose().rename(columns={0 : "count"})
    df["percentage"] = (df["count"] / nb) * 100
    df = df.reset_index(drop=False).rename(columns={"index" : "distance"})
    g = bar(data=df,x="distance",y="percentage",color=color,ax=ax)
    for idx,p in enumerate (g.patches) : 
        g.annotate(str(df.iloc[idx]["count"]),(p.get_x() + p.get_width() / 2 , p.get_height()))
    g.set_xlabel("Distance to TSS (kb)")
    g.set_ylabel("Genomic region (%)")
    g.set_title("Binned by absolute distance to TSS")

def graph_absolute_dist_tss(test:str or pd.DataFrame,regdom:str or pd.DataFrame,ax=None,color="#325fa8") -> None : 
    """
    This function allows the creation of a barplot of the absolute
    distance between the peaks and the TSS of the associated gene(s). 

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame 
        Regulatory domain of all genes in the genome 
    ax : 
        Define the position of the plot in a figure

    Returns
    -------
    None

    Exemples 
    --------
    Example available here: https://github.com/theislab/greatpy/blob/main/notebooks/07_plot.ipynb    
   
    """
    res = {"0:5": [0],"5:50": [0],"50:500": [0],">500": [0]}
    nb = 0
    dist = gp.tl.get_dist_to_tss(test,regdom)

    for i in dist.values() : 
        for j in i : 
            j = abs(j)
            if j < 5000 :
                res["0:5"][0] += 1
            elif j < 50000 :
                res["5:50"][0] += 1
            elif j < 500000 :
                res["50:500"][0] += 1
            else :
                res[">500"][0] += 1
            nb += 1
    df = pd.DataFrame(res).transpose().rename(columns={0:"count"})
    df["percentage"] = (df["count"]/nb)*100
    df = df.reset_index(drop=False).rename(columns={"index":"distance"})
    g = bar(data=df,x="distance",y="percentage",color=color,ax=ax)
    for idx,p in enumerate (g.patches) : 
        g.annotate(str(df.iloc[idx]["count"]),(p.get_x() + p.get_width() / 2 , p.get_height()))
    g.set_xlabel("Absolute distance to TSS (kb)")
    g.set_ylabel("Genomic region (%)")
    g.set_title("Binned by absolute distance to TSS")
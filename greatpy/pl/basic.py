import pandas as pd 
from numpy import log
import greatpy as gp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams

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
        sb.scatterplot(data=great_df,x = f"-log({colname_x})",y = f"-log({colname_y})",ax = ax).set_title(title)
    else : 
        sb.scatterplot(data=great_df,x=colname_x,y=colname_y,ax=ax).set_title(title)


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

    g = sb.barplot(data = nb,x="number",y="percentage",ax=ax,color=color)
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
    g = sb.barplot(data=df,x="distance",y="percentage",color=color,ax=ax)
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
    g = sb.barplot(data=df,x="distance",y="percentage",color=color,ax=ax)
    for idx,p in enumerate (g.patches) : 
        g.annotate(str(df.iloc[idx]["count"]),(p.get_x() + p.get_width() / 2 , p.get_height()))
    g.set_xlabel("Absolute distance to TSS (kb)")
    g.set_ylabel("Genomic region (%)")
    g.set_title("Binned by absolute distance to TSS")


def scale_data_5_75(data):
    mind = np.min(data)
    maxd = np.max(data)
    
    if maxd == mind:
        maxd=maxd+1
        mind=mind-1
        
    drange = maxd - mind
    return ((((data - mind)/drange*0.70)+0.05)*100)


def plot_enrich(data, n_terms=20, color="cool", save=False):
    """
    This function allows the creation of a dotplot of the enrichment
    GO term in the inputs datas  

    Parameters
    ----------
    data : pd.DataFrame
        Results for greatpy 
    n_terms : int 
        the number of term to be shown 
    save : bool 
        Is the plot should be save 

    Returns
    -------
    None

    Exemples 
    --------
   
    """
    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError('Please input a Pandas Dataframe output by gprofiler.')
        
    if not np.all([term in data.columns for term in ['p_value', 'name', 'intersection_size']]):
        raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))
    
    data_to_plot = data.iloc[:n_terms,:].copy()
    data_to_plot['go.id'] = data_to_plot.index

    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap=color, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'bold'})

    sb.set(style="whitegrid")

    path = plt.scatter(x='recall', y="name", c='p_value', cmap=color, 
                    norm=colors.LogNorm(min_pval, max_pval), 
                    data=data_to_plot, linewidth=1, edgecolor="grey", 
                    s=[(i+10)**1.5 for i in data_to_plot['scaled.overlap']])
    ax = plt.gca()
    ax.invert_yaxis()

    ax.set_ylabel('')
    ax.set_xlabel('Gene ratio', fontsize=14, fontweight='bold')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    #Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick)/6).astype(int)
    
    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick-1
    
    ticks_vals = [10**i for i in range(max_tick, min_tick-1, -tick_step)]
    ticks_labs = ['$10^{'+str(i)+'}$' for i in range(max_tick, min_tick-1, -tick_step)]

    #Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([0.8, 0.15, 0.03, 0.4])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0,0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs)
    cbar.set_label("Adjusted p-value", fontsize=14, fontweight='bold')

    #Size legend
    min_olap = data_to_plot['intersection_size'].min()
    max_olap = data_to_plot['intersection_size'].max()
    olap_range = max_olap - min_olap
    
    #Note: approximate scaled 5, 25, 50, 75 values are calculated
    #      and then rounded to nearest number divisible by 5
    size_leg_vals = [np.ceil(i/5)*5 for i in 
                        [min_olap, min_olap+(20/70)*olap_range, min_olap+(45/70)*olap_range, max_olap]]
    size_leg_scaled_vals = scale_data_5_75(size_leg_vals)

    
    l1 = plt.scatter([],[], s=(size_leg_scaled_vals[0]+10)**1.5, edgecolors='none', color='black')
    l2 = plt.scatter([],[], s=(size_leg_scaled_vals[1]+10)**1.5, edgecolors='none', color='black')
    l3 = plt.scatter([],[], s=(size_leg_scaled_vals[2]+10)**1.5, edgecolors='none', color='black')
    l4 = plt.scatter([],[], s=(size_leg_scaled_vals[3]+10)**1.5, edgecolors='none', color='black')

    labels = [str(int(i)) for i in size_leg_vals]

    leg = plt.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False, fontsize=12,
                    handlelength=1, loc = 'center left', borderpad = 1, labelspacing = 1.4,
                    handletextpad=2, title='Gene overlap', scatterpoints = 1,  bbox_to_anchor=(-2, 1.5), 
                    facecolor='black')

    if save:
        plt.savefig(save, dpi=300, format='pdf', bbox_inches='tight')

    plt.show()
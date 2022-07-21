import pandas as pd 
import numpy as np 
from math import lgamma, log, exp,fabs,inf
from scipy.special import comb
import greatpy as gp

def get_association(test,regdom): 
    """
    Function allowing from a file of genomic regions from CHIPseq 
    and a file of genomic regulatory domains to determine the names 
    of genes associated with at least one genomic region 

    Parameters
    ----------
    test : pd.dataFrame
        df of the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : list
        list of gene associated with at least with one test peak
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> get_association(test,regdom)        
        ['RNF223']
    
    """
    res = []
    for i in range(test.shape[0]) :
        currTest = test.iloc[i]
        regdom_curr_test = regdom.loc[(regdom["Chr"] == currTest["Chr"])].sort_values("Chr_Start")
        regdom_curr_test = regdom_curr_test.loc[
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"])) | # regdom overlap totally test 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"])) | # test overlap totally regdom 
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_Start"])) | # regdom overlap not totally test on left side 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"]) & (regdom_curr_test["Chr_Start"] <= currTest["Chr_End"])) # regdom overlap not totally test on right side 
            ] 
        res = res + list(regdom_curr_test["Name"])
    return list(dict.fromkeys(res))

def len_regdom(regdom:pd.DataFrame): 
    """
    Function to calculate for each gene name from regdom the
     size of the regulatory region for this gene in the genome 

    Parameters
    ----------    
    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    dict
        dictionary in which each key corresponds to a gene name 
        from regdom and the value is the size of the regulatory 
        region for that gene
        
    Exemples 
    --------
    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']}))

    >>> len_regdom(regdom)
        {'RNF223': 80097, 'C1orf159': 52710}

    """
    test = regdom["Chr_End"]-regdom["Chr_Start"]
    return pd.DataFrame({"len":list(test)},index = regdom["Name"]).to_dict()["len"]

def number_of_hit(test,regdom): 
    """ 
    Function to calculate the number of hits from several 
    genomic regions and the file describing the regulatory regions

    Parameters
    ----------
    test : pd.dataFrame
        df of the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    nb : int
        number of hit 
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> number_of_hit(test,regdom)        
        1
    
    """
    nb = 0
    regdom = regdom[["Chr","Chr_Start","Chr_End"]]
    for i in range(test.shape[0]) : 
        chrom = test.iat[i,0]
        start = test.iat[i,1]
        end = test.iat[i,2]
        regdom_np = regdom["Chr"].to_numpy()
        reg_start = regdom["Chr_Start"].to_numpy()
        reg_end = regdom["Chr_End"].to_numpy()
        Chr_reduce = np.where(regdom_np == chrom)
        reg_start = np.take(reg_start,Chr_reduce,axis = 0)[0]
        reg_end = np.take(reg_end,Chr_reduce,axis = 0)[0]

        if any((reg_start <= start) & (reg_end >= end)):  
            nb += 1
    return nb

def betacf(a,b,x): 
    """ Used by betai: Evaluates continued fraction for incomplete beta function """
    maxit = 10000
    eps = 3.0e-7 
    fpmin = 1.0e-30
    qab = a+b
    qap = a+1
    qam = a-1
    c = 1
    d = 1-qab*x/qap
    if fabs(d) < fpmin :
        d = fpmin
    d = 1/d 
    h = d
    for m in range(1,maxit+1) : 
        m2 = 2*m
        aa = m*(b-m)*x / ((qam+m2) * (a+m2))
        d = 1.0+aa*d
        if (fabs(d) < fpmin) : 
            d = fpmin
        c = 1.0+aa/c
        if (fabs(c) < fpmin):
            c = fpmin
        d = 1.0/d
        h *= d*c
        aa = -(a+m) * (qab+m)*x / ((a+m2) * (qap+m2))
        d = 1.0+aa*d  
        if (fabs(d) < fpmin):
            d = fpmin
        c = 1.0+aa/c
        if (fabs(c) < fpmin):
            c = fpmin
        d = 1.0/d
        dell = d*c
        h *= dell
        if (fabs(dell-1.0) < eps):
            break
    if (m > maxit):
        print("a or b too big, or MAXIT too small in betacf")
        return False
    return h

def betai(a,b,x):
    """Returns the incomplete beta function Ix(a, b)."""
    if x < 0 or x > 1 : 
        # print("bad x in routine betai")
        return False
    if x == 0 or x == 1 : 
        bt = 0.0
    else : 
        bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x))
    if x < (a+1) / (a+b+2) : 
        return bt * betacf(a,b,x)/a
    return 1 - bt*betacf(b,a,1-x)/b

def get_binom_pval(n:int,k:int,p:float) -> float:
    """
    This function allows to calculate the binomial probability 
    of obtaining k in a set of size n and whose probability is p 

    Parameters
    ----------
    n : int
        Number of genomic region in the test set 
    k : int 
        Number of test genomic regions in the regulatory domain of a gene with annotation
    p : float
        Percentage of genome annotated

    Returns
    -------
    float
        binomial probability
        
    Exemples 
    --------
    >>> get_binom_pval(100,2,0.2)
        0.9999999947037065
    
    """
    if k == 0 : return 1
    else : return betai(k,n-k+1,p)

def hypergeom_pmf(N, K, n, k):
    """
    Function to calculate the probability mass function for hypergeometric distribution

    Parameters
    ----------
    N : int
        Total number of gene in the genome
    K : int 
        Number of genes in the genome with annotation
    n : int
        Number of gene in the test set
    k : int
        Number of genes in the test gene set with annotation

    Returns
    -------
    float
        proability mass function
        
    Exemples 
    --------
    >>> hypergeom_pmf(100,10,30,1)
        0.11270773995748315
    
    """
    Achoosex = comb(K,k,exact=True) 
    NAchoosenx = comb(N-K, n-k,exact=True) 
    Nchoosen = comb(N,n,exact=True) 
    return ((Achoosex)*NAchoosenx)/Nchoosen 

def hypergeom_cdf(N, K, n, k):
    """
    Function to calculate the cumulative density funtion for hypergeometric distribution

    Parameters
    ----------
    N : int
        Total number of gene in the genome
    K : int 
        Number of genes in the genome with annotation
    n : int
        Number of gene in the test set
    k : int
        Number of genes in the test gene set with annotation

    Returns
    -------
    float
        Cumulative density function
        
    Exemples 
    --------
    >>> hypergeom_cdf(100,10,30,1)
        0.9770827595419788
    
    """
    return np.sum([hypergeom_pmf(N, K, n, x) for x in range(k,min(K,n)+1)])


def get_nb_asso_per_region(test,regdom) : 
    """
    Function allowing from a file of genomic regions from CHIPseq 
    and a file of genomic regulatory domains to determine number of peaks 
    associated with each gene in the regulatory domain. 

    Parameters
    ----------
    test : str
        path of the file with the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : str
        path of the file with the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the number of associated genes per genomic region : key = associated gene, value = number of peaks associated with the gene 
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> get_association(test,regdom)        
        {'RNF223':2}
    
    """
    res = {}
    test,regdom,_,_ = gp.tl.GREAT.loader(test,regdom,None,None) 
    for i in range(test.shape[0]) :
        currTest = test.iloc[i]
        regdom_curr_test = regdom.loc[(regdom["Chr"] == currTest["Chr"])].sort_values("Chr_Start")
        regdom_curr_test = regdom_curr_test.loc[
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"])) | # regdom overlap totally test 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"])) | # test overlap totally regdom 
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_Start"])) | # regdom overlap not totally test on left side 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"]) & (regdom_curr_test["Chr_Start"] <= currTest["Chr_End"])) # regdom overlap not totally test on right side 
            ] 
        res[i] = regdom_curr_test.shape[0]
    return res

def get_dist_to_tss(test,regdom) : 
    """
    Function allowing from a file of genomic regions from CHIPseq 
    and a file of genomic regulatory domains to determine the distance from peaks 
    to the transcription start site of the associated gene

    Parameters
    ----------
    test : str
        path of the file with the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : str
        path of the file with the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the distance from tss to the associated genes : key = number of the input, value = distance from peaks to tss of associated genes 
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> get_association(test,regdom)        
        {'RNF223':[-22278]}
    
    """
    res = {}
    test,regdom,_,_ = gp.tl.GREAT.loader(test,regdom,None,None) 
    for i in range(test.shape[0]) :
        currTest = test.iloc[i]
        mean_pos_test = (currTest["Chr_End"] + currTest["Chr_Start"])/2
        regdom_curr_test = regdom.loc[(regdom["Chr"] == currTest["Chr"])].sort_values("Chr_Start")
        regdom_curr_test = regdom_curr_test.loc[
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"])) | # regdom overlap totally test 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"])) | # test overlap totally regdom 
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_Start"])) | # regdom overlap not totally test on left side 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"]) & (regdom_curr_test["Chr_Start"] <= currTest["Chr_End"])) # regdom overlap not totally test on right side 
            ] 
        res[i] = []
        for j in range(regdom_curr_test.shape[0]) : 
            res[i].append(regdom_curr_test.iloc[j]["tss"] - mean_pos_test)
    return res
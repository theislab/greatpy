import pandas as pd
from math import lgamma, log, exp,fabs,inf
from scipy.stats import hypergeom
from scipy.special import comb
from statsmodels.stats.multitest import multipletests
import dask.dataframe as dd 
import numpy as np 

pd.options.display.float_format = '{:12.5e}'.format

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
    for i in range(test.shape[0]):
        currTest = test.iloc[i]
        regdom_curr_test = regdom.loc[regdom["Chr"] == currTest["Chr"]].sort_values("Chr_Start")
        regdom_inf = regdom_curr_test.loc[regdom_curr_test["tss"] <= currTest["Chr_Start"]]
        regdom_sup = regdom_curr_test.loc[regdom_curr_test["tss"] >= currTest["Chr_End"]]
        try : 
            if regdom_inf.iloc[-1]["Name"] not in res : 
                res.append (regdom_inf.iloc[-1]["Name"])
        except :
            pass
        try :
            if regdom_sup.iloc[0]["Name"] not in res :
                res.append(regdom_sup.iloc[0]["Name"])
        except : 
            pass
    return res

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
    return pd.DataFrame({"len":list(test)},index=regdom["Name"]).to_dict()["len"]

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
    nb=0
    regdom = regdom[["Chr","Chr_Start","Chr_End"]]
    for i in range(test.shape[0]): 
        chrom = test.iat[i,0]
        start = test.iat[i,1]
        end = test.iat[i,2]
        regdom_np = regdom["Chr"].to_numpy()
        reg_start = regdom["Chr_Start"].to_numpy()
        reg_end = regdom["Chr_End"].to_numpy()
        Chr_reduce = np.where(regdom_np == chrom)
        reg_start = np.take(reg_start,Chr_reduce,axis=0)[0]
        reg_end = np.take(reg_end,Chr_reduce,axis=0)[0]

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
    for m in range(1,maxit+1): 
        m2 = 2*m
        aa = m*(b-m)*x/((qam+m2)*(a+m2))
        d = 1.0+aa*d
        if (fabs(d) < fpmin) : 
            d = fpmin
        c = 1.0+aa/c
        if (fabs(c) < fpmin):
            c = fpmin
        d = 1.0/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
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
        print("bad x in routine betai")
        return False
    if x == 0 or x == 1 : 
        bt = 0.0
    else : 
        bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x))
    if x < (a+1)/(a+b+2) : 
        return bt*betacf(a,b,x)/a
    return 1-bt*betacf(b,a,1-x)/b

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
    Achoosex = comb(K,k) if comb(K,k) != inf else 1e-308
    NAchoosenx = comb(N-K, n-k) if comb(N-K, n-k) != inf else 1e-308
    Nchoosen = comb(N,n) if comb(N,n) != inf else 1e-308
    return ((Achoosex)*NAchoosenx)/Nchoosen if Nchoosen > 1e-308 and (Achoosex)*NAchoosenx != 0.0 else 1e-308

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

class GREAT: 
    def loader(test_data:str or pd.DataFrame,regdom_file:str or pd.DataFrame,chr_size_file:str or pd.DataFrame,annotation_file:str or pd.DataFrame):
        """
        This function is used to load all datasets needed for the enrichment calculation

        Parameters
        ----------
        test_data : str or pd.DataFrame
            Genomic set of peaks to be tested
        regdom_file : str or pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size_file : str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : str or pd.DataFrame
            Table with the annotation of each gene in the genome

        Returns
        -------
        test_data : pd.DataFrame
            Genomic set of peaks to be tested in the good format 
        regdom : pd.DataFrame
            Regulatory domain of all genes in the genome in the good format
        size : pd.DataFrame
            Table with the size of each chromosome in the good format
        ann : pd.DataFrame
            Table with the annotation of each gene in the genome in the good format
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader(
        ...    "../../data/human/test_genomic_region.bed",
        ...    "../../data/human/regulatory_domain.bed",
        ...    "../../data/human/chr_size.bed",
        ...    "../../data/human/ontologies.csv"
        ...    )

        >>> test.head()
            |    | Chr   |   Chr_Start |   Chr_End |
            |---:|:------|------------:|----------:|
            |  0 | chr1  |     1052028 |   1052049 |
            |  1 | chr1  |     1065512 |   1065533 |
            |  2 | chr1  |     1067375 |   1067397 |
            |  3 | chr1  |     1068083 |   1068119 |
            |  4 | chr1  |    10520283 |  10520490 |

        >>> regdom.head()
            |    | Chr   |   Chr_Start |   Chr_End | Name      |   tss | Strand   |
            |---:|:------|------------:|----------:|:----------|------:|:---------|
            |  0 | chr1  |           0 |     22436 | MIR6859-1 | 17436 | -        |
            |  1 | chr1  |       16436 |     22436 | MIR6859-2 | 17436 | -        |
            |  2 | chr1  |       16436 |     22436 | MIR6859-3 | 17436 | -        |
            |  3 | chr1  |       16436 |     28370 | MIR6859-4 | 17436 | -        |
            |  4 | chr1  |       22436 |     34370 | WASH7P    | 29370 | -        |

        >>> size.head()
            |    | Chrom   |      Size |
            |---:|:--------|----------:|
            |  0 | chr1    | 248956422 |
            |  1 | chr2    | 242193529 |
            |  2 | chr3    | 198295559 |
            |  3 | chr4    | 190214555 |
            |  4 | chr5    | 181538259 |
            
        >>> ann.head()
            |    | id         | name                                                   | symbol        |
            |---:|:-----------|:-------------------------------------------------------|:--------------|
            |  0 | GO:0003924 | GTPase activity                                        | DNAJC25-GNG10 |
            |  1 | GO:0007186 | G protein-coupled receptor signaling pathway           | DNAJC25-GNG10 |
            |  2 | GO:0003723 | RNA binding                                            | NUDT4B        |
            |  3 | GO:0005829 | cytosol                                                | NUDT4B        |
            |  4 | GO:0008486 | diphosphoinositol-polyphosphate diphosphatase activity | NUDT4B        |
        
        """

        if type(regdom_file) == str:
            regdom = pd.read_csv(regdom_file,sep="\t",comment="#",
                        names=["Chr", "Chr_Start", "Chr_End","Name","tss","Strand"],dtype={"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64","Name":"object","tss":"int64","Strand":"object"})
        else:
            regdomfile = regdom_file.iloc[:,:6]
            colname = list(regdomfile.columns)
            try : 
                regdom_file = regdom_file.rename(columns={colname[0]:"Chr",colname[1]:"Chr_Start",colname[2]:"Chr_End","Name":"Name","tss":"tss","Strand":"Strand"})
            except :
                print("Error in the format of the regdom file")
                print("The regdom file must have the following columns : Chr, Chr_Start, Chr_End, Name, tss, Strand")
                return False 

        if type(test_data) == str : 
            test_data = pd.read_csv(test_data,sep="\t",comment="#",usecols=[0,1,2],
                            names=["Chr", "Chr_Start", "Chr_End"],dtype={"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64"})
        else : 
            test_data = test_data.iloc[:,:3]
            colname = list(test_data.columns)
            try : 
                test_data = test_data.rename(columns={colname[0]:"Chr",colname[1]:"Chr_Start",colname[2]:"Chr_End"})
            except : 
                print("Error in test dataframe, please check your input")
                print("Columns should be : chr...(type object), start(type int), end(type int)")
                return False
        
        if type(chr_size_file) == str :
            size = pd.read_csv(chr_size_file,sep="\t",comment="#",
                            names=["Chrom","Size"],dtype={"Chrom":"object", "Size":"int64"})
        else :
            chr_size_file = chr_size_file.iloc[:,:2]
            colname = list(chr_size_file.columns)
            try : 
                chr_size_file = chr_size_file.rename(columns={colname[0]:"Chrom",colname[1]:"Size"})
            except : 
                print("Error in the format of the chr_size file")
                print("The chr_size file must have the following columns : Chrom, Size")
                return False
        if type(annotation_file) == str : 
            dask_df = dd.read_csv(annotation_file,sep=";",  comment = "#",
                            dtype={"ensembl":"object","id":"object","name":"object","ontology.group":"object","gene.name":"object","symbol":"object"},
                            usecols=["id","name","symbol"],low_memory=False)
            ann = dask_df.compute()
            ann = ann[ann['id'].str.match('^GO.*')== True]
        else : 
            ann = annotation_file.iloc[:,:4]
            colname = list(ann.columns)
            try : 
                ann = ann.rename(columns={colname[0]:"id",colname[1]:"name",colname[3]:"symbol"})
            except : 
                print("Error in the format of the annotation file")
                print("The annotation file must have the following columns : id, name, symbol")
                return False
        return test_data,regdom,size,ann

    def __enrichment_binom_and_hypergeom(test,regdom,size,ann,asso) : 
        """
        This private function is used to compute the enrichment of the test data using the binomial test and the hypergeometric test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list 
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the binomial test and the the hypergeometric test
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.____enrichment_binom_and_hypergeom(
        ...    test = test,
        ...    regdom = regdom,
        ...    size = size,
        ...    ann = ann,
        ...    asso = get_association(test,regdom)
        ...    )

        >>> enrichment.head()
            |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |
            |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|
            | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  |
            | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  |
            | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  |
            | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 |
            | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  |
        """
        # Init Great
        res = {}
        hit={}

        # init Hypergeom
        hypergeom_gene_set=len(asso) # get the number of genes in the test gene set.
        hypergeom_total_number_gene = regdom.shape[0] #get the number of genes in the genome.

        # Init binom 
        n_binom = test.shape[0]# get the number of genomic region in the test set
        total_nu = size["Size"].sum()# get the total number of nucleotides in the genome

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = len_regdom(regdom)# get the length of each regulatory domain 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id'].isin([i])]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation

                if i not in list(hit.keys()) : 
                    hit[i] = number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom=hit[i]
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom))
            res.update({elem[2]:[ elem[3],get_binom_pval(n_binom,elem[0],elem[1]/total_nu), hypergeom_cdf(hypergeom_total_number_gene,elem[4],hypergeom_gene_set,elem[5]) ] for elem in tmp})
        return pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"binom_p_value",2:"hypergeom_p_value"}).replace(0,np.nan).sort_values(by="binom_p_value")
    
    def __enrichment_binom(test,regdom,size,ann,asso):
        """
        This private function is used to compute the enrichment of the test data using the binomial test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list 
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the binomial test
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.____enrichment_binom(
        ...    test = test,
        ...    regdom = regdom,
        ...    size = size,
        ...    ann = ann,
        ...    asso = get_association(test,regdom)
        ...    )

        >>> enrichment.head()
            |            | go_term                                                          |   binom_p_value |
            |:-----------|:-----------------------------------------------------------------|----------------:|
            | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |
            | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |
            | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |
            | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |
            | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |
        """ 
        # Init Great
        res = {}
        hit = {}

        # Init binom 
        n_binom = test.shape[0]# get the number of genomic region in the test set
        total_nu = size["Size"].sum()# get the total number of nucleotides in the genome
        
        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = len_regdom(regdom)# get the length of each regulatory domain 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            tmp=[]
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id'].isin([i])]
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]

                if i not in list(hit.keys()) : 
                    hit[i] = number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom=hit[i]
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"]))
            res.update({elem[2]:[ elem[3],get_binom_pval(n_binom,elem[0],elem[1]/total_nu) ] for elem in tmp})
        return pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"binom_p_value"}).sort_values(by="binom_p_value").sort_values(by="binom_p_value")

    def __enrichment_hypergeom(test,regdom,ann,asso): 
        """
        This private function is used to compute the enrichment of the test data using the hypergeometric test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list 
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the hypergeometric test
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.____enrichment_hypergeom(
        ...    test = test,
        ...    regdom = regdom,
        ...    ann = ann,
        ...    asso = get_association(test,regdom)
        ...    )

        >>> enrichment.head()   
        """ 
        # Init Great
        res = {}

        # Init hypergeom
        hypergeom_total_number_gene = regdom.shape[0] #get the number of genes in the genome
        hypergeom_gene_set = len(asso) # get the number of genes in the test gene set.

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"] == name]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id']==i]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation                
                tmp.append((i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom)) 
            res.update({elem[0]:[ elem[1], hypergeom_cdf(hypergeom_total_number_gene,elem[2],hypergeom_gene_set,elem[3]) ] for elem in tmp}) 
        return pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"hypergeom_p_value"}).replace(0,3e-308).sort_values(by="hypergeom_p_value")


    def enrichment(test_file,regdom_file,chr_size_file, annotation_file, binom=True,hypergeom=True):
        """
        This function is a wrapper of the 3 private methods: 
        * GREAT.__enrichment_binom_and_hypergeom 
        * GREAT.__enrichment_binom 
        * GREAT.__enrichment_hypergeom

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list 
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the hypergeometric test
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.enrichment(
        ...    test = test,
        ...    regdom = regdom,
        ...    ann = ann,
        ...    asso = get_association(test,regdom),
        ...    binom=True,
        ...    hypergeom=True
        ...    )
        >>> enrichment.head()
            |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |
            |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|
            | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  |
            | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  |
            | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  |
            | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 |
            | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  |

        >>> enrichment = GREAT.enrichment(
        ...    test = test,
        ...    regdom = regdom,
        ...    ann = ann,
        ...    asso = get_association(test,regdom),
        ...    binom=True,
        ...    hypergeom=False
        ...    )
        >>> enrichment.head()
            |            | go_term                                                          |   binom_p_value |
            |:-----------|:-----------------------------------------------------------------|----------------:|
            | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |
            | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |
            | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |
            | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |
            | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |

        >>> enrichment = GREAT.enrichment(
        ...    test = test,
        ...    regdom = regdom,
        ...    ann = ann,
        ...    asso = get_association(test,regdom),
        ...    binom=False,
        ...    hypergeom=True
        ...    )
        >>> enrichment.head()

        """
        if not binom and not hypergeom : 
            return False
        
        test,regdom,size,ann = GREAT.loader(test_file,regdom_file,chr_size_file, annotation_file)
        asso = get_association(test,regdom)# get the name of the regulatory domain associated to each genomic region in the test set


        if binom and hypergeom : 
            return GREAT.__enrichment_binom_and_hypergeom(test,regdom,size,ann,asso)

        elif binom : 
            return GREAT.__enrichment_binom(test,regdom,size,ann,asso)

        else : 
              return GREAT.__enrichment_hypergeom(test,regdom,ann,asso)

    def set_bonferroni(self,alpha:float=0.05): 
        """
        This function create new columns in the dataframe with the Bonferroni correction

        Parameters
        ----------
        alpha : float
            alpha value for the Bonferroni correction
        
        Returns
        -------
        pd.DataFrame
            dataframe new columns with the Bonferroni correction for each p-value
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.enrichment(test = test,regdom = regdom,ann = ann,asso = get_association(test,regdom),binom=True,hypergeom=True)
        >>> bonferroni = GREAT.set_bonferroni(enrichment,alpha=0.05)
        >>> bonferroni.head()

        """
        for col in self.columns: 
            if col in ["binom_p_value","hypergeom_p_value"] : 
                col_split = col.split("_")
                self[f"{col_split[0]}_bonferroni"] = multipletests(self[col], alpha=alpha, method='bonferroni')[1]
        return self 

    def set_fdr(self,alpha:float=0.05) : 
        """
        This function create new columns in the dataframe with the fdr correction

        Parameters
        ----------
        alpha : float
            alpha value for the fdr correction
        
        Returns
        -------
        pd.DataFrame
            dataframe new columns with the fdr correction for each p-value
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.enrichment(test = test,regdom = regdom,ann = ann,asso = get_association(test,regdom),binom=True,hypergeom=True)
        >>> fdr = GREAT.set_fdr(enrichment,alpha=0.05)
        >>> fdr.head()

        """
        for col in self.columns: 
            if col in ["binom_p_value","hypergeom_p_value"] : 
                col_split = col.split("_")
                # self[f"{col_split[0]}_fdr"] = fdrcorrection(self[col], alpha=alpha)[1]
                self[f"{col_split[0]}_fdr"] = multipletests(self[col], alpha=alpha, method='fdr_bh')[1]
        return self 

    def set_threshold(self,colname:str, alpha:int=0.05) : 
        """
        This function allows to delete rows according to the p-value of the column taken as argument. By default the alpha value is 0.05

        Parameters
        ----------
        alpha : float
            alpha value for the fdr correction
        
        Returns
        -------
        pd.DataFrame
            dataframe with the rows deleted according to the p-value threshold
            
        Exemples 
        --------
        >>> test,regdom,size,ann = GREAT.loader("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv")
        >>> enrichment = GREAT.enrichment(test = test,regdom = regdom,ann = ann,asso = get_association(test,regdom),binom=True,hypergeom=True)
        >>> significant = GREAT.set_threshold(enrichment,colname="binom_p_value",alpha=0.05)
        >>> significant.head()

        """
        if colname in self.columns: 
            self = self.loc[self[colname]<=alpha]
        return self 

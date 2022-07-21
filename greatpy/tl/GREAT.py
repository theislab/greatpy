import pandas as pd
from statsmodels.stats.multitest import multipletests
import dask.dataframe as dd 
import numpy as np 
import greatpy as gp 
pd.options.display.float_format = '{:12.5e}'.format

class GREAT: 
    def loader(test_data:None or str or pd.DataFrame,regdom_file:None or str or pd.DataFrame,chr_size_file:None or str or pd.DataFrame,annotation_file:None or str or pd.DataFrame):
        """
        This function is used to load all datasets needed for the enrichment calculation

        Parameters
        ----------
        test_data : None or str or pd.DataFrame
            Genomic set of peaks to be tested
        regdom_file : None or str or pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size_file : None or str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : None or str or pd.DataFrame
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
        ...    "../../data/tests/test_data/input/02_srf_hg38.bed",
        ...    "../../data/human/hg38/regulatory_domain.bed",
        ...    "../../data/human/hg38/chr_size.bed",
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

        if type(regdom_file) == str :
            regdom = pd.read_csv(regdom_file,sep = "\t",comment="#",
                        names = ["Chr", "Chr_Start", "Chr_End","Name","tss","Strand"],
                        dtype = {"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64","Name":"object","tss":"int64","Strand":"object"})
        elif type(regdom_file) == pd.DataFrame :
            regdom = regdom_file.iloc[:,:6]
            colname = list(regdom.columns)
            try : 
                regdom = regdom.rename(columns = {colname[0]:"Chr",colname[1]:"Chr_Start",colname[2]:"Chr_End",colname[3]:"Name",colname[4]:"tss",colname[5]:"Strand"})
            except :
                print("Error in the format of the regdom file")
                print("The regdom file must have the following columns : Chr, Chr_Start, Chr_End, Name, tss, Strand")
                return False 
        else : 
            regdom = regdom_file

        if type(test_data) == str : 
            test_data = pd.read_csv(test_data,sep = "\t",comment = "#",usecols = [0,1,2],
                            names = ["Chr", "Chr_Start", "Chr_End"],
                            dtype = {"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64"})
        elif type(test_data) == pd.DataFrame : 
            test_data = test_data.iloc[:,:3]
            colname = list(test_data.columns)
            try : 
                test_data = test_data.rename(columns={colname[0]:"Chr",colname[1]:"Chr_Start",colname[2]:"Chr_End"})
            except : 
                print("Error in test dataframe, please check your input")
                print("Columns should be : chr...(type object), start(type int), end(type int)")
                return False
        else :
            pass

        if type(chr_size_file) == str :
            size = pd.read_csv(chr_size_file,sep = "\t",comment = "#",
                            names = ["Chrom","Size"],
                            dtype = {"Chrom":"object", "Size":"int64"})
        elif type(chr_size_file) == pd.DataFrame :
            size = chr_size_file.iloc[:,:2]
            colname = list(size.columns)
            try : 
                size = size.rename(columns = {colname[0]:"Chrom",colname[1]:"Size"})
            except : 
                print("Error in the format of the chr_size file")
                print("The chr_size file must have the following columns : Chrom, Size")
                return False
        else :
            size = chr_size_file

        if type(annotation_file) == str : 
            dask_df = dd.read_csv(annotation_file,sep = ";",  comment = "#",
                            dtype = {"ensembl":"object","id":"object","name":"object","ontology.group":"object","gene.name":"object","symbol":"object"},
                            usecols = ["id","name","symbol"],low_memory = False)
            ann = dask_df.compute()
            ann = ann[ann['id'].str.match('^GO.*')== True]
        elif type(annotation_file) == pd.DataFrame : 
            ann = annotation_file.iloc[:,:4]
            colname = list(ann.columns)
            try : 
                ann = ann.rename(columns={colname[0]:"id",colname[1]:"name",colname[3]:"symbol"})
            except : 
                print("Error in the format of the annotation file")
                print("The annotation file must have the following columns : id, name, symbol")
                return False
        else :
            ann = annotation_file

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
        >>> test,regdom,size,ann = GREAT.loader("../../data/tests/test_data/input/02_srf_hg38.bed","../../data/human/hg38/regulatory_domain.bed","../../data/human/hg38/chr_size.bed","../../data/human/ontologies.csv")
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
        hit = {}

        # init Hypergeom
        hypergeom_gene_set = len(asso) # get the number of genes in the test gene set.
        hypergeom_total_number_gene = regdom.shape[0] #get the number of genes in the genome.

        # Init binom 
        n_binom = test.shape[0]# get the number of genomic region in the test set
        total_nu = size["Size"].sum()# get the total number of nucleotides in the genome

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = gp.tl.len_regdom(regdom)# get the length of each regulatory domain 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())) : 
                gene_imply = ann[ann['id'].isin([i])]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation

                if i not in list(hit.keys()) : 
                    hit[i] = gp.tl.number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom = hit[i]
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom))
            res.update({elem[2]:[ elem[3],gp.tl.get_binom_pval(n_binom,elem[0],elem[1]/total_nu),elem[0]/(elem[1]/total_nu), gp.tl.hypergeom_cdf(hypergeom_total_number_gene,elem[4],hypergeom_gene_set,elem[5]),(elem[5]*hypergeom_total_number_gene)/(hypergeom_gene_set*elem[4]) ] for elem in tmp})
        return pd.DataFrame(res).transpose().rename(columns = {0:"go_term",1:"binom_p_value",2:"binom_fold_enrichment",3:"hypergeom_p_value",4:"hypergeometric_fold_enrichment"}).replace(0,np.nan).dropna().sort_values(by = "binom_p_value")
    
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
        >>> test,regdom,size,ann = GREAT.loader("../../data/tests/test_data/input/02_srf_hg38.bed","../../data/human/hg38/regulatory_domain.bed","../../data/human/hg38/chr_size.bed","../../data/human/ontologies.csv")
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
        len_on_chr = gp.tl.len_regdom(regdom)# get the length of each regulatory domain 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())) : 
                gene_imply = ann[ann['id'].isin([i])]
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]

                if i not in list(hit.keys()) : 
                    hit[i] = gp.tl.number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom = hit[i]
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"]))
            res.update({elem[2]:[ elem[3],gp.tl.get_binom_pval(n_binom,elem[0],elem[1]/total_nu),elem[0]/(elem[1]/total_nu) ] for elem in tmp})
        return pd.DataFrame(res).transpose().rename(columns = {0:"go_term",1:"binom_p_value",2:"binom_fold_enrichment"}).sort_values(by = "binom_p_value")

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
            for i in (list(id.unique())) : 
                gene_imply = ann[ann['id']==i]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation                
                tmp.append((i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom)) 
            res.update({elem[0]:[ elem[1], gp.tl.hypergeom_cdf(hypergeom_total_number_gene,elem[2],hypergeom_gene_set,elem[3]),(elem[3]*hypergeom_total_number_gene)/(hypergeom_gene_set*elem[2]) ] for elem in tmp}) 
        return pd.DataFrame(res).transpose().rename(columns = {0:"go_term",1:"hypergeom_p_value",2:"hypergeometric_fold_enrichment"}).replace(0,np.nan).dropna().sort_values(by = "hypergeom_p_value")


    def enrichment(test_file: str or pd.DataFrame,regdom_file: str or pd.DataFrame,chr_size_file: str or pd.DataFrame, annotation_file: str or pd.DataFrame, binom=True,hypergeom=True):
        """
        This function is a wrapper of the 3 private methods: 
        * GREAT.__enrichment_binom_and_hypergeom 
        * GREAT.__enrichment_binom 
        * GREAT.__enrichment_hypergeom

        Parameters
        ----------
        test_file : str or pd.DataFrame
            Genomic set of peaks to be tested
        regdom_file : str or pd.DataFrame 
            Regulatory domain of all genes in the genome 
        chr_size_file : str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : str or pd.DataFrame
            Table with the annotation of each gene in the genome
        binom : bool (default True)
            If True, the binomial test is used.
        hypergeom : bool (default True)
            If True, the hypergeometric test is used.

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
        ...    chr_size_file = size,
        ...    ann = ann,
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
        asso = gp.tl.get_association(test,regdom)# get the name of the regulatory domain associated to each genomic region in the test set


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
        >>> enrichment = great.tl.GREAT.enrichment("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv",binom=True,hypergeom=True)
        >>> bonferroni = GREAT.set_bonferroni(enrichment,alpha=0.05)
        >>> bonferroni.head()
        ...    |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |   binom_bonferroni |   hypergeom_bonferroni |
        ...    |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|-------------------:|-----------------------:|
        ...    | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  |        3.0754e-10  |                      1 |
        ...    | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  |        2.87384e-07 |                      1 |
        ...    | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  |        2.87384e-07 |                      1 |
        ...    | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 |        3.74892e-07 |                      1 |
        ...    | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  |        7.68931e-07 |                      1 |

        """
        for col in self.columns : 
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
        >>> enrichment = great.tl.GREAT.enrichment("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv",binom=True,hypergeom=True)
        >>> fdr = GREAT.set_fdr(enrichment,alpha=0.05)
        >>> fdr.head()
        ...    |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |   binom_fdr |   hypergeom_fdr |
        ...    |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|------------:|----------------:|
        ...    | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  | 3.0754e-10  |       0.0913909 |
        ...    | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  | 9.3723e-08  |       0.0913909 |
        ...    | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  | 9.3723e-08  |       0.0913909 |
        ...    | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 | 9.3723e-08  |       0.0913909 |
        ...    | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  | 1.53786e-07 |       0.0913909 |

        """
        for col in self.columns : 
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
        >>> enrichment = great.tl.GREAT.enrichment("../../data/human/test_genomic_region.bed", "../../data/human/regulatory_domain.bed", "../../data/human/chr_size.bed", "../../data/human/ontologies.csv",binom=True,hypergeom=True)
        >>> enrichment.shape[0]
        ...    594

        >>> significant = GREAT.set_threshold(enrichment,colname="binom_p_value",alpha=0.05)
        >>> significant.shape[0]
        ...    310

        """
        if colname in self.columns: 
            self = self.loc[self[colname]<=alpha]
        return self 
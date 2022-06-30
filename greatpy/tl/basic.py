from anndata import AnnData
import pandas as pd
from math import lgamma, log, exp,fabs
pd.options.display.float_format = '{:12.5e}'.format
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests,fdrcorrection
from scipy.stats import hypergeom as hg 


def basic_tool(adata: AnnData) -> int:
    """Run a tool on the AnnData object."""
    print("Implement a tool to run on the AnnData object.")
    return 0


max_extension = 1000000  
basal_upstream = 5000
basal_downstream = 1000

def file_reader(path:str,type_f:str):
    """
    path of the file and type of data in the file 

    type should be : 
        tss if : "Chr","tss","strand","name"
        chr_size if : "Chr","Chr_Size" 
    """
    if type_f not in ["tss","chr_size"]: 
        print("type should be tss or chr_size")
        return False

    elif type_f == "tss": 
        data=pd.read_csv(path,sep="\t",comment="#",names=["Chr","tss","Strand","name"])
        return data

    elif type_f == "chr_size": 
        data=pd.read_csv(path,sep="\t",comment="#",names=["Chr","Size"])
        return data

def validate_input(association): 
    if association != "OneCloset" and association != "TwoCloset" and association != "Basalplusextention" : 
        print("Association rule should be OneCloset or TwoCloset Basalplusextention")
        return False
    if max_extension < 0: 
        print(f"Maximum extension must be a non-negative integer: {max_extension}")
        return False
    if basal_upstream < 0 : 
        print(f"Basal upstream must be a non-negative integer: {basal_upstream}")
        return False
    if basal_downstream < 0 : 
        print(f"Basal downstream must be a non-negative integer: {basal_downstream}")
        return False
    return True

def write_Regdom(regdom:pd.DataFrame,file_name):
    f = open(file_name,"w")
    f.write("#chr\tChrStart\tChrEnd\tname\ttss\tstrand\n")
    for i in range(regdom.shape[0]): 
        curr = regdom.iloc[i]
        chr = curr["Chr"];start = curr["Chr_Start"];end = curr["Chr_End"];name = curr["name"];tss = curr["tss"];strand = curr["Strand"]
        f.write(f"{chr}\t{start}\t{end}\t{name}\t{tss}\t{strand}\n")
    f.close()

def create_basal_plus_extension_regdom(regdom:pd.DataFrame,maximumExtension,basalUp,basalDown,chr_size):
    prev = curr = next = 0
    chr_strat_end = []
    start = []
    end = []
    for i in range(regdom.shape[0]):
        curr = regdom.iloc[i]
        chr = curr["Chr"]
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        tmp = curr["tss"]
        if curr["Strand"] == "+": 
            curr_chr_start = max(0,tmp-basalUp)
            curr_chr_end = min(curr_chr_size,tmp+basalDown)
        elif curr["Strand"] == "-": 
            curr_chr_start = max(0,tmp-basalDown)
            curr_chr_end = min(curr_chr_size,tmp+basalUp)
        elif curr["Strand"] == "." : 
            print("Invalid_Input : Impossible to create a basal expression regdom if you have not specify the strand")
            return False
        else : 
            err = curr["Strand"]
            print (f"Invalid input : strand should be '+' or '-'. Line {i} : strand = {err}")
            return False
        chr_strat_end.append([curr_chr_start,curr_chr_end])

    for i in range(regdom.shape[0]):
        curr = regdom.iloc[i]
        chr = curr["Chr"]
        
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        if i != regdom.shape[0]-1 : 
            next = regdom.iloc[i+1]
        else : next = 0


        tmp_start = max(0,curr["tss"]-maximumExtension)
        basal_start = chr_strat_end[i][0]
        tmp_start = min(basal_start,tmp_start)
        if type(prev) != int and prev["Chr"] == curr["Chr"]:
            if prev["Strand"] == "+": 
                prev_end = prev["tss"]+basalDown
            else : 
                prev_end = prev["tss"]+basalUp
            tmp_start = min(basal_start,max(prev_end,tmp_start))


        tmp_end = min(curr_chr_size,curr["tss"]+max_extension)
        basal_end = chr_strat_end[i][1]
        
        tmp_end = max(basal_end,tmp_end)
        if type(next) != int and next["Chr"] == curr["Chr"]:
            if next["Strand"] == "+": 
                nextStart = next["tss"]-basalUp
            else : 
                nextStart = next["tss"]-basalDown
            tmp_end = max(basal_end,min(nextStart,tmp_end))
            
        prev = regdom.iloc[i]
        start.append(int(tmp_start))
        end.append(int(tmp_end))
    regdom["Chr_Start"] = start
    regdom["Chr_End"] = end
    return regdom

def create_Two_Closet_Regdom(regdom,max_extension,chr_size):
    return create_basal_plus_extension_regdom(regdom,max_extension,0,0,chr_size)

def create_one_closet_regdom(regdom:pd.DataFrame,maximum_extension,chr_size:pd.DataFrame):
    prev = curr = next = 0
    start = []
    end = []

    for i in range(regdom.shape[0]):
        curr = regdom.iloc[i]
        chr = curr["Chr"]
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        if i < regdom.shape[0]-1 : next = regdom.iloc[i+1]
        else : next = 0

        if (type(prev) == int and prev == 0) or (type(prev) != int and prev["Chr"] != curr["Chr"]) : prev = 0 
        if (type(next) == int and next == 0) or (type(next) != int and next["Chr"] != curr["Chr"]): next = 0

        tmp_start = max(0,curr["tss"]-maximum_extension)
        if type(prev) != int : 
            middle = (curr["tss"]+prev["tss"])//2
            tmp_start = max(tmp_start,middle)
        
        tmp_end = min(curr["tss"]+maximum_extension,curr_chr_size)
        if type(next) != int  :
            middle = (curr["tss"]+next["tss"])//2
            tmp_end = min(tmp_end,middle)
        
        start.append(tmp_start)
        end.append(tmp_end)
        prev = curr
        
    regdom["Chr_Start"] = start
    regdom["Chr_End"] = end
    return regdom

def create_regdom(tss_file,chr_sizes_file,association_rule,out_path): 
    if not validate_input(association_rule): 
        print("Invalid input")
        return False
    df = file_reader(tss_file,'tss')
    
    df = df.sort_values(["Chr","tss","Strand","name"])

    chr_size = file_reader(chr_sizes_file,"chr_size") 

    if association_rule == "OneCloset" : 
        out = create_one_closet_regdom(df,max_extension,chr_size)
    elif association_rule == "TwoCloset" : 
        out = create_Two_Closet_Regdom(df,max_extension,chr_size)
    elif association_rule == "Basalplusextention" : 
        out = create_basal_plus_extension_regdom(df,max_extension,basal_upstream,basal_downstream,chr_size)
    else : 
        return False
    out = out.astype({"Chr_Start":int,"Chr_End":int})
    out = out.reindex(["Chr","Chr_Start","Chr_End","name","tss","Strand"],axis=1)

    write_Regdom(out,out_path) 
    return out

def get_association(test,regdom): 
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

def len_regdom(regdom): 
    test = regdom["Chr_End"]-regdom["Chr_Start"]
    return pd.DataFrame({"len":list(test)},index=regdom["Name"]).to_dict()["len"]

def number_of_hit(test,regdom): 
    nb=0
    for i in range(test.shape[0]): 
        curr_test = test.iloc[i]
        regdom_reduce = regdom.loc[regdom["Chr"] == curr_test["Chr"]]
        if regdom_reduce[(regdom_reduce["Chr_Start"] <= curr_test["Chr_Start"]) & (regdom_reduce["Chr_End"] >= curr_test["Chr_End"])].shape[0] > 0 : 
            nb += 1
    return nb

def betacf(a,b,x): 
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

def get_binom_pval(n,k,p):
    if k == 0 : return 1
    else : return betai(k,n-k+1,p)

def enrichment(test:str or pd.DataFrame,regdom_file,chr_size_file,annotation,binom=True,hypergeom=True,alpha=0.05,correction=("fdr",0.05),sort_by=None): 
    # Data import 
    if not binom and not hypergeom : 
        return False
    
    regdom = pd.read_csv(regdom_file,sep="\t",comment="#",
                    names=["Chr", "Chr_Start", "Chr_End","Name","tss","Strand"],dtype={"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64","Name":"object","tss":"int64","Strand":"object"})

    if type(test) == str : 
        test = pd.read_csv(test,sep="\t",comment="#",
                        names=["Chr", "Chr_Start", "Chr_End"],dtype={"Chr":"object", "Chr_Start":"int64", "Chr_End":"int64"})
    else : 
        test = test.iloc[:,:3]
        colname = list(test.columns)
        try : 
            test = test.rename(columns={colname[0]:"Chr",colname[1]:"Chr_Start",colname[2]:"Chr_End"})
        except : 
            print("Error in test dataframe, please check your input")
            print("Columns should be : chr...(type object), start(type int), end(type int)")
            return False
    
    size = pd.read_csv(chr_size_file,sep="\t",comment="#",
                    names=["Chrom","Size"],dtype={"Chrom":"object", "Size":"int64"})

    ann = pd.read_csv(annotation,sep=";",  
                    names=["ensembl","id","name","ontology.group","gene.name","symbol"],dtype={"ensembl":"object","id":"object","name":"object","ontology.group":"object","gene.name":"object","symbol":"object"},
                    usecols=["id","name","gene.name","symbol"],low_memory=True)
    ann = ann[ann['id'].str.match('^GO.*')== True]

    if binom and hypergeom : 
        # Init Great
        res = {}
        hypergeom_total_number_gene = regdom.shape[0] #get the number of genes in the genome.
        n_binom = test.shape[0]# get the number of genomic region in the test set
        total_nu = size["Size"].sum()# get the total number of nucleotides in the genome
        asso = get_association(test,regdom)# get the name of the regulatory domain associated to each genomic region in the test set
        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = len_regdom(regdom)# get the length of each regulatory domain 

        # init Hypergeom
        hypergeom_gene_set=len(asso) # get the number of genes in the test gene set.

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"]==name]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id']==i]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["gene.name"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation
                k_binom = number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom))
            res.update({elem[2]:[ elem[3],get_binom_pval(n_binom,elem[0],elem[1]/total_nu), sum([hg.pmf(i,hypergeom_total_number_gene,hypergeom_gene_set,elem[4]) for i in range(elem[5],min(elem[4],hypergeom_gene_set)+1)]) ] for elem in tmp})
        
        df= pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"binom_p_value",2:"hypergeom_p_value"})
        if correction == (0,0) or correction[0] not in ['bonferroni','fdr'] or correction[1] >= 1 or correction[1] <= 0: 
            return df.sort_values(by=sort_by) if sort_by != None else df 

        elif correction[0] == "bonferroni" : 
            df["binom_bonferroni_correction"] = multipletests(df["binom_p_value"], alpha=correction[1], method='bonferroni')[1]
            df["hypergeom_bonferroni_correction"] = multipletests(df["hypergeom_p_value"], alpha=correction[1], method='bonferroni')[1]
            df = df.loc[df["binom_bonferroni_correction"] <= alpha]

        elif correction[0] == "fdr" : 
            df["binom_fdr_correction"] = fdrcorrection(df["binom_p_value"], alpha=correction[1])[1]
            df["hypergeom_fdr_correction"] = fdrcorrection(df["hypergeom_p_value"], alpha=correction[1])[1]
            df = df.loc[df["binom_fdr_correction"] <= alpha]
        return df.sort_values(by=sort_by) if sort_by != None else df 

    elif binom : 
        # Init Great
        res = {}
        n_binom = test.shape[0]# get the number of genomic region in the test set
        total_nu = size["Size"].sum()# get the total number of nucleotides in the genome
        asso = get_association(test,regdom)# get the name of the regulatory domain associated to each genomic region in the test set
        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = len_regdom(regdom)# get the length of each regulatory domain 

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"] == name]
            id = ann_name_gene["id"]
            tmp=[]
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id']==i]
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["gene.name"]))]
                k_binom = number_of_hit(test,curr_regdom)# get the number of test genomic regions in the regulatory domain of a gene with annotation
                nb_binom = sum([len_on_chr[i] for i in curr_regdom["Name"]])# get the portion of the genome in the regulatory domain of a gene with annotation
                
                tmp.append((k_binom,nb_binom,i,gene_imply.iloc[0]["name"]))
            res.update({elem[2]:[ elem[3],get_binom_pval(n_binom,elem[0],elem[1]/total_nu) ] for elem in tmp})
        df= pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"binom_p_value"}).sort_values(by="binom_p_value")
        if correction == (0,0) or correction[0] not in ['bonferroni','fdr'] or correction[1] >= 1 or correction[1]<=0: 
            return df.sort_values(by=sort_by) if sort_by != None else df 

        elif correction[0] == "bonferroni" : 
            df["binom_bonferroni_correction"] = multipletests(df["binom_p_value"], alpha=correction[1], method='bonferroni')[1]
            df = df.loc[df["binom_bonferroni_correction"]<=alpha]

        elif correction[0] == "fdr" :
            df["binom_fdr_correction"] = fdrcorrection(df["binom_p_value"], alpha=correction[1])[1]
            df = df.loc[df["binom_fdr_correction"]<=alpha] 
        return df.sort_values(by=sort_by) if sort_by != None else df 

    else : 
        # Init Great
        res = {}
        hypergeom_total_number_gene = regdom.shape[0] #get the number of genes in the genome.
        asso = get_association(test,regdom)# get the name of the regulatory domain associated to each genomic region in the test set
        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]#reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association 
        len_on_chr = len_regdom(regdom)# get the length of each regulatory domain 

        # init Hypergeom
        hypergeom_gene_set = len(asso) # get the number of genes in the test gene set.

        #Compute for all associating gene and for each GO id associated with the gene the probability. 
        for name in asso :
            ann_name_gene = ann[ann["symbol"] == name]
            id = ann_name_gene["id"]
            tmp = []
            for i in (list(id.unique())): 
                gene_imply = ann[ann['id']==i]
                K_hypergeom = gene_imply.shape[0] # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["gene.name"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[0] # get the number of genes in the test gene set with annotation                
                tmp.append((i,gene_imply.iloc[0]["name"],K_hypergeom,k_hypergeom)) 
            res.update({elem[0]:[ elem[1], sum([hg.pmf(i,hypergeom_total_number_gene,hypergeom_gene_set,elem[2]) for i in range(elem[3],min(elem[2],hypergeom_gene_set)+1)]) ] for elem in tmp}) 

        df = pd.DataFrame(res).transpose().rename(columns={0:"go_term",1:"hypergeom_p_value"}).sort_values(by="hypergeom_p_value")
        
        if correction == (0,0) or correction[0] not in ['bonferroni','fdr'] or correction[1] >= 1 or correction[1]<=0: 
            df = df.loc[df["hypergeom_p_value"] <= alpha]
            return df.sort_values(by=sort_by) if sort_by != None else df  

        elif correction[0] == "bonferroni" : 
            df["hypergeom_bonferroni_correction"] = multipletests(df["hypergeom_p_value"], alpha=correction[1], method='bonferroni')[1]
            df = df.loc[df["hypergeom_bonferroni_correction"] <= alpha]

        elif correction[0] == "fdr" : 
            df["hypergeom_fdr_correction"] = fdrcorrection(df["hypergeom_p_value"], alpha=correction[1])[1]
            df = df.loc[df["hypergeom_fdr_correction"] <= alpha]

        return df.sort_values(by=sort_by) if sort_by != None else df  
from anndata import AnnData
import pandas as pd
from math import lgamma, log, exp,fabs
pd.options.display.float_format = '{:12.5e}'.format


def basic_tool(adata: AnnData) -> int:
    """Run a tool on the AnnData object."""
    print("Implement a tool to run on the AnnData object.")
    return 0


maxExtension = 1000000  
basalUpstream = 5000
basalDownstream = 1000

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

def validate_input(Association): 
    if Association!="OneCloset" and Association!="TwoCloset" and Association!="Basalplusextention" : 
        print("Association rule should be OneCloset or TwoCloset Basalplusextention")
        return False
    if maxExtension < 0: 
        print(f"Maximum extension must be a non-negative integer: {maxExtension}")
        return False
    if basalUpstream < 0 : 
        print(f"Basal upstream must be a non-negative integer: {basalUpstream}")
        return False
    if (basalDownstream < 0) : 
        print(f"Basal downstream must be a non-negative integer: {basalDownstream}")
        return False
    return True

def write_Regdom(regdom:pd.DataFrame,file_name):
    f=open(file_name,"w")
    f.write("#chr\tChrStart\tChrEnd\tname\ttss\tstrand\n")
    for i in range(regdom.shape[0]): 
        curr=regdom.iloc[i]
        chr=curr["Chr"];start=curr["Chr_Start"];end=curr["Chr_End"];name=curr["name"];tss=curr["tss"];strand=curr["Strand"]
        f.write(f"{chr}\t{start}\t{end}\t{name}\t{tss}\t{strand}\n")
    f.close()

def create_Basal_Plus_Extension_Regdom(regdom:pd.DataFrame,maximumExtension,basalUp,basalDown,Chr_size):
    prev=curr=next=0
    chr_strat_end=[]
    start=[]
    end=[]
    for i in range(regdom.shape[0]):
        curr=regdom.iloc[i]
        chr=curr["Chr"]
        curr_chr_size=int(Chr_size.loc[Chr_size["Chr"]==chr,"Size"])
        tmp=curr["tss"]
        if curr["Strand"]=="+": 
            curr_chr_start=max(0,tmp-basalUp)
            curr_chr_end=min(curr_chr_size,tmp+basalDown)
        elif curr["Strand"]=="-": 
            curr_chr_start=max(0,tmp-basalDown)
            curr_chr_end=min(curr_chr_size,tmp+basalUp)
        elif curr["Strand"]=="." : 
            print("Invalid_Input : Impossible to create a basal expression regdom if you have not specify the strand")
            return False
        else : 
            err=curr["Strand"]
            print (f"Invalid input : strand should be '+' or '-'. Line {i} : strand = {err}")
            return False
        chr_strat_end.append([curr_chr_start,curr_chr_end])

    for i in range(regdom.shape[0]):
        curr=regdom.iloc[i]
        chr=curr["Chr"]
        
        curr_chr_size=int(Chr_size.loc[Chr_size["Chr"]==chr,"Size"])
        if i != regdom.shape[0]-1 : 
            next=regdom.iloc[i+1]
        else : next=0


        tmpStrat=max(0,curr["tss"]-maximumExtension)
        basalStart=chr_strat_end[i][0]
        tmpStrat=min(basalStart,tmpStrat)
        if type(prev)!=int and prev["Chr"]==curr["Chr"]:
            if prev["Strand"]=="+": 
                prevEnd=prev["tss"]+basalDown
            else : 
                prevEnd=prev["tss"]+basalUp
            tmpStrat=min(basalStart,max(prevEnd,tmpStrat))


        tmpEnd=min(curr_chr_size,curr["tss"]+maxExtension)
        basalEnd=chr_strat_end[i][1]
        
        tmpEnd=max(basalEnd,tmpEnd)
        if type(next)!=int and next["Chr"]==curr["Chr"]:
            if next["Strand"]=="+": 
                nextStart=next["tss"]-basalUp
            else : 
                nextStart=next["tss"]-basalDown
            tmpEnd=max(basalEnd,min(nextStart,tmpEnd))
            
        prev=regdom.iloc[i]
        start.append(int(tmpStrat))
        end.append(int(tmpEnd))
    regdom["Chr_Start"]=start
    regdom["Chr_End"]=end
    return regdom

def create_Two_Closet_Regdom(regdom,maxExtension,Chr_size):
    return create_Basal_Plus_Extension_Regdom(regdom,maxExtension,0,0,Chr_size)

def create_One_Closet_Regdom(regdom:pd.DataFrame,maximumExtension,Chr_size:pd.DataFrame):
    prev=curr=next=0
    start=[]
    end=[]

    for i in range(regdom.shape[0]):
        curr=regdom.iloc[i]
        chr=curr["Chr"]
        curr_chr_size=int(Chr_size.loc[Chr_size["Chr"]==chr,"Size"])
        if i < regdom.shape[0]-1 : next=regdom.iloc[i+1]
        else : next=0

        if (type(prev)==int and prev==0) or (type(prev)!=int and prev["Chr"] != curr["Chr"]) : prev=0 
        if (type(next)==int and next==0) or (type(next)!=int and next["Chr"] != curr["Chr"]): next=0

        tmp_start=max(0,curr["tss"]-maximumExtension)
        if type(prev)!=int : 
            middle = (curr["tss"]+prev["tss"])//2
            tmp_start=max(tmp_start,middle)
        
        
        tmp_end=min(curr["tss"]+maximumExtension,curr_chr_size)
        if type(next)!=int  :
            middle = (curr["tss"]+next["tss"])//2
            tmp_end=min(tmp_end,middle)
        
        start.append(tmp_start)
        end.append(tmp_end)
        prev=curr
        
    regdom["Chr_Start"]=start
    regdom["Chr_End"]=end
    return regdom

def create_Regdom(tssFn,chromSizesFn,AssociationRule,outFn): 
    if not validate_input(AssociationRule): 
        print("Invalid input")
        return False
    df = file_reader(tssFn,'tss')
    
    df= df.sort_values(["Chr","tss","Strand","name"])

    chr_size=file_reader(chromSizesFn,"chr_size") 

    if AssociationRule == "OneCloset" : 
        out=create_One_Closet_Regdom(df,maxExtension,chr_size)
    elif AssociationRule == "TwoCloset" : 
        out=create_Two_Closet_Regdom(df,maxExtension,chr_size)
    elif AssociationRule == "Basalplusextention" : 
        out=create_Basal_Plus_Extension_Regdom(df,maxExtension,basalUpstream,basalDownstream,chr_size)
    else : 
        return False
    out=out.astype({"Chr_Start":int,"Chr_End":int})
    out= out.reindex(["Chr","Chr_Start","Chr_End","name","tss","Strand"],axis=1)

    write_Regdom(out,outFn) 
    return out

def get_association(test,regdom): 
    res = []
    for i in range(test.shape[0]):
        currTest=test.iloc[i]
        regdom_curr_test = regdom.loc[regdom["Chr"]==currTest["Chr"]].sort_values("Chr_Start")
        regdom_inf = regdom_curr_test.loc[regdom_curr_test["tss"]<=currTest["Chr_Start"]]
        regdom_sup = regdom_curr_test.loc[regdom_curr_test["tss"]>=currTest["Chr_End"]]
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
    res={}
    for i in range(regdom.shape[0]): 
        res[regdom.iloc[i]["Name"]]=regdom.iloc[i]["Chr_End"]-regdom.iloc[i]["Chr_Start"]
    return res

def hit(test,regdom): 
    nb=0
    for i in range(test.shape[0]): 
        currTest=test.iloc[i]
        regdom=regdom.loc[regdom["Chr"]==currTest["Chr"]]
        for j in range(regdom.shape[0]): 
            currRegdom=regdom.iloc[j]
            if currTest["Chr_Start"] >= currRegdom["Chr_Start"] and currTest["Chr_End"] <= currRegdom["Chr_End"] :
                nb+=1
    return nb

def betacf(a,b,x): 
    MAXIT = 10000
    EPS = 3.0e-7 
    FPMIN = 1.0e-30
    qab=a+b
    qap=a+1
    qam=a-1
    c=1
    d=1-qab*x/qap
    if fabs(d) < FPMIN :
        d = FPMIN
    d=1/d 
    h=d
    for m in range(1,MAXIT+1): 
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0+aa*d
        if (fabs(d) < FPMIN) : 
            d=FPMIN
        c=1.0+aa/c
        if (fabs(c) < FPMIN):
            c=FPMIN
        d=1.0/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.0+aa*d  
        if (fabs(d) < FPMIN):
            d=FPMIN
        c=1.0+aa/c
        if (fabs(c) < FPMIN):
            c=FPMIN
        d=1.0/d
        dell=d*c
        h *= dell
        if (fabs(dell-1.0) < EPS):
            break
    if (m > MAXIT):
        print("a or b too big, or MAXIT too small in betacf")
        return False
    return h

def betai(a,b,x):
    if x < 0 or x > 1 : 
        print("bad x in routine betai")
        return False
    if x == 0 or x == 1 : 
        bt=0.0
    else : 
        bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x))
    if x < (a+1)/(a+b+2) : 
        return bt*betacf(a,b,x)/a
    return 1-bt*betacf(b,a,1-x)/b

def get_Binom_Pval(n,k,p):
    if k == 0 : return 1
    else : return betai(k,n-k+1,p)

def calculBinomP(test,regdomFn,Chr_sizeFn,annotation): 
    regdom=pd.read_csv(regdomFn,sep="\t",comment="#",names=["Chr", "Chr_Start", "Chr_End","Name","tss","Strand"])
    test=pd.read_csv(test,sep="\t",comment="#",names=["Chr", "Chr_Start", "Chr_End"])
    n = test.shape[0] # get the number of genomic region in the test set
    size=pd.read_csv(Chr_sizeFn,sep="\t",comment="#",names=["Chrom","Size"])
    G = size["Size"].sum() # get the total number of nucleotides in the genome
    ann = pd.read_csv(annotation,sep="\t",names=["ensembl","id","name","ontology.group","gene.name","symbol"],low_memory=False)
    ann = ann[ann['id'].str.match('^GO.*')== True]

    res={}
    asso= get_association(test,regdom) # get the name of the regulatory domain associated to each genomic region in the test set
    
    ann_red = ann[ann["symbol"].isin(asso)]
    regdom = regdom[regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))]
    len_on_chr=len_regdom(regdom) # get the length of each regulatory domain
    for name in asso :
        ann_name_gene = ann[ann["symbol"]==name]
        id=ann_name_gene["id"]
        tmp=[]
        for i in (list(id.unique())): 
            gene_imply=ann[ann['id']==i]
            curr_regdom=regdom.loc[regdom["Name"].isin(list(gene_imply["gene.name"]))]
            k = hit(test,curr_regdom) # get the number of test genomic regions in the regulatory domain of a gene with annotation
            nb=sum([len_on_chr[i] for i in curr_regdom["Name"]]) # get the portion of the genome in the regulatory domain of a gene with annotation
            tmp.append((k,nb,i,gene_imply.iloc[0]["name"]))
        
        res.update({elem[2]:[elem[3],get_Binom_Pval(n,elem[0],elem[1]/G)] for elem in tmp})# if get_Binom_Pval(n,elem[0],elem[1]/G)<0.05 and get_Binom_Pval(n,elem[0],elem[1]/G)!=0})

    return pd.DataFrame(res).transpose().rename(columns={0:"GO_term",1:"P-value"}).sort_values(by="P-value")
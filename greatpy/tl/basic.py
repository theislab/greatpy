from anndata import AnnData
import pandas as pd
from math import lgamma, log, exp


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

def get_range_tree_of_regdom(regdom): 
    return regdom[["Chr","Chr_Start","Chr_End"]]

def get_Total_Non_Gap_Bases(antigap):
    retval=0
    for i in range(antigap.shape[0]): 
        retval+= antigap.iloc[i]["Chr_End"]-antigap.iloc[i]["Chr_Start"]
    return retval

def genomeOverlapSize(ranges,chr,start,end):
    """ 
    Function to get the number of intersection between the range and an element which we know the chr number, the start and stop position on the chromosome
    """
    # ranges=ranges.loc[ranges["Chr"]==chr & ranges["Chr_Start"]>start & ranges["Chr_End"]<end] # don't work I don't understand
    ranges=ranges.loc[ranges["Chr"]==chr]
    ranges=ranges.loc[ranges["Chr_Start"]>start]
    ranges=ranges.loc[ranges["Chr_End"]<end]
    for i in range(ranges.shape[0]): 
        if ranges.iloc[i]["Chr_Start"] > start : 
            start = ranges.iloc[i]["Chr_Start"]
        elif ranges.iloc[i]["Chr_End"] < end : 
            end = ranges.iloc[i]["Chr_End"]
    return end-start

def get_anotated_Non_Gap_Bases(ranges,antigap):
    """Function to get all of the intersection of the antigap and the ranges data """
    retval=0
    for i in range(antigap.shape[0]):
        currAntigap=antigap.iloc[i]
        retval+=genomeOverlapSize(ranges,currAntigap["Chr"],currAntigap["Chr_Start"],currAntigap["Chr_End"])
    return retval

def betacf(a,b,x): 
    MAXIT = 10000
    EPS = 3.0e-7 
    FPMIN = 1.0e-30
    qab=a+b
    qap=a+1
    qam=a-1
    c=1
    d=1-qab*x/qap
    if abs(d) < FPMIN :
        d = FPMIN
    d=1/d 
    h=d
    for m in range(1,MAXIT+1): 
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0+aa*d
        if (abs(d) < FPMIN) : 
            d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN):
            c=FPMIN
        d=1.0/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.0+aa*d  
        if (abs(d) < FPMIN):
            d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN):
            c=FPMIN
        d=1.0/d
        dell=d*c
        h *= dell
        if (abs(dell-1.0) < EPS):
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
        return bt*betacf(a,b,x)
    return 1-bt*betacf(b,a,1-x)/b

def get_Binom_Pval(n,k,p):
    if k == 0 : return 1
    else : return betai(k,n-k+1,p)

def calculBinomP(regdomFn,antigapFn,totalRegions,hitRegions) : 
    df=pd.read_csv(regdomFn,sep="\t",comment="#",names=["Chr", "Chr_Start", "Chr_End","Name","tss","Strand"])
    ranges=get_range_tree_of_regdom(df)
    antigap=pd.read_csv(antigapFn,sep="\t",comment="#",names=["Chr", "Chr_Start", "Chr_End","Characteristics"])
    
    total_Non_Gap_Bases=get_Total_Non_Gap_Bases(antigap)
    anotated_Non_Gap_Bases=get_anotated_Non_Gap_Bases(ranges,antigap)

    annotation_Weight=anotated_Non_Gap_Bases/total_Non_Gap_Bases

    binomP = get_Binom_Pval(totalRegions,hitRegions,annotation_Weight)
    return binomP 
from anndata import AnnData
import pandas as pd
from math import lgamma, log, exp,fabs,inf
pd.options.display.float_format = '{:12.5e}'.format
from scipy.stats import hypergeom
from scipy.special import comb
from statsmodels.stats.multitest import multipletests,fdrcorrection
from scipy.stats import hypergeom as hg 
import dask.dataframe as dd 
import cython
import numpy as np 


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

    if out_path != None : 
        write_Regdom(out,out_path) 
    return out


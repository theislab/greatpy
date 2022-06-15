from anndata import AnnData
import pandas as pd


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
        print("type should be tsv,tss,BED or chr_size")
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